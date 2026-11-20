"""Functions shared by the analyses of library pools.

A pool is sequenced over a serial dilution. Diluting it should divide the viral counts by
the dilution factor while leaving the neutralization standard alone, so the odds of a read
coming from the neutralization standard rise along a line of slope `IDEAL_SLOPE` against
the dilution factor on log-log axes for as long as that holds. The composition of the pool
can be read only over the dilutions where it does: a well too concentrated for its counts
to fall with dilution measures the assay rather than the pool.

These functions locate that linear range, choose the wells to read, and turn the counts of
those wells into each unit's share of the pool.

"""

import math
from dataclasses import dataclass

import pandas

# What diluting the pool should do to the neutralization-standard odds, the standard's
# counts over the viral counts: halving the pool doubles the odds. The odds are what the
# report plots, so the fit is on them and a slope reads with the same sign as the line drawn.
# This is what the assay should do rather than a judgement call, so it is not configurable.
IDEAL_SLOPE = 1.0

# Where the pool starts responding to dilution is judged over windows of this many consecutive
# dilutions, and a window's fitted slope counts as ideal when it is within this tolerance of
# `IDEAL_SLOPE`. The windows are short and local on purpose: a window running to the most
# dilute end of the series averages the wells that behave with those that do not.
DEFAULT_LINEARITY_WINDOW_DILUTIONS = 3
DEFAULT_LINEARITY_SLOPE_TOL = 0.3

# How many standard errors of its own slope widen the tolerance for a window. Widening is only
# ever consulted on the explicit-wells path: a window admitted by it is reported as such but
# never chooses the wells, since the latitude grows with the noise and so admits a
# saturated, noisy window more readily than a cleaner one.
DEFAULT_LINEARITY_SLOPE_Z = 3

# the unit a pool can be balanced by, and so the unit an equal share is an equal share of
POOL_BY_OPTIONS = ("strain", "barcode", "stock")


@dataclass(frozen=True)
class AnalysisWells:
    """The wells a pool's composition is read from, and where its titer comes from.

    Attributes:
        `wells` (list)
            The ``(plate, well)`` of each well to read, all at one dilution.
        `dilution` (int or float)
            The dilution those wells are at.
        `titer_window` (list or None)
            The dilutions of the window the titer is computed over, or `None` where no
            strict window gives one and the titer is therefore not reported.
        `strict` (bool)
            Whether that window was within tolerance without being widened. False where
            the tolerance had to be widened to place the range at all, and where there is
            no window.
        `wells_given` (bool)
            Whether the wells were configured rather than located from the fit.

    """

    wells: list
    dilution: float
    titer_window: list
    strict: bool
    wells_given: bool


def log_log_slope(x, y):
    """Least-squares fit of log10(`y`) on log10(`x`).

    Returns the slope, its standard error, and the intercept, which is what places the
    fitted line where a plot draws it.

    """
    # named here rather than left to `math.log10`, which raises without saying which of the
    # two it was given or that a logarithm was what needed the value
    for name, values in [("x", x), ("y", y)]:
        nonpositive = sorted(v for v in values if not v > 0)
        if nonpositive:
            raise ValueError(
                f"cannot fit on a logarithmic scale: {name} holds {nonpositive}, and a "
                "logarithm is defined only for a positive value"
            )
    x = [math.log10(v) for v in x]
    y = [math.log10(v) for v in y]
    if len(x) < 3:
        raise ValueError("need at least three points to fit a slope and its error")
    mean_x, mean_y = sum(x) / len(x), sum(y) / len(y)
    sxx = sum((v - mean_x) ** 2 for v in x)
    if sxx == 0:
        raise ValueError("cannot fit a slope to a single distinct dilution")
    slope = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y, strict=True)) / sxx
    sse = sum(
        (b - mean_y - slope * (a - mean_x)) ** 2 for a, b in zip(x, y, strict=True)
    )
    return slope, math.sqrt(sse / (len(x) - 2) / sxx), mean_y - slope * mean_x


def fit_dilution_windows(
    wells,
    window_dilutions=DEFAULT_LINEARITY_WINDOW_DILUTIONS,
    slope_tol=DEFAULT_LINEARITY_SLOPE_TOL,
    slope_z=DEFAULT_LINEARITY_SLOPE_Z,
):
    """Fit the slope of every window of `window_dilutions` consecutive dilutions.

    Args:
        `wells` (pandas.DataFrame)
            The wells to fit, with ``dilution_factor`` and ``neut_standard_odds``
            columns. Only wells that passed QC should be passed, as a well that did not is
            not trusted to measure the pool.
        `window_dilutions` (int)
            How many consecutive dilutions a window spans. At least three, since fewer cannot
            give a slope and its error, and odd, so that a window has a middle dilution.
        `slope_tol` (float)
            How far a window's slope may be from `IDEAL_SLOPE` and still count as ideal.
        `slope_z` (float)
            Standard errors of the slope by which that tolerance is widened.

    Returns:
        List of dicts, one per window from the most concentrated end, each with the window's
        ``window`` of dilutions, its ``slope``, ``slope_se`` and ``intercept``, which
        together place its fitted line, its ``deviation`` from
        `IDEAL_SLOPE`, and whether it is ``within_tol`` and ``within_tol_widened``. Empty
        where the series is shorter than one window.

    """
    if window_dilutions < 3 or window_dilutions % 2 == 0:
        raise ValueError(
            f"'linearity_window_dilutions' must be an odd number of at least 3, not "
            f"{window_dilutions}: fewer than three dilutions cannot give a slope and its "
            "error, and an even window has no middle dilution to read the composition from"
        )
    # `tolist` rather than `unique`, so the dilutions are Python numbers: a `numpy` one
    # carries its type into the reported window and writes as `np.int64(128)`
    dilutions = sorted(set(wells["dilution_factor"].tolist()))
    windows = []
    for i in range(len(dilutions) - window_dilutions + 1):
        window = dilutions[i : i + window_dilutions]
        in_window = wells[wells["dilution_factor"].isin(window)]
        slope, slope_se, intercept = log_log_slope(
            in_window["dilution_factor"], in_window["neut_standard_odds"]
        )
        deviation = abs(slope - IDEAL_SLOPE)
        windows.append(
            {
                "window": window,
                "slope": slope,
                "slope_se": slope_se,
                # with the slope, where this window's own fitted line is drawn
                "intercept": intercept,
                "deviation": deviation,
                "within_tol": deviation <= slope_tol,
                "within_tol_widened": deviation <= slope_tol + slope_z * slope_se,
            }
        )
    return windows


def analysis_wells(windows, wells, given=None):
    """Choose the wells the composition is read from, and the window the titer comes from.

    Without `given`, the wells are the middle dilution of the first strictly linear window,
    rather than its most concentrated edge, so that there is a dilution of margin on either
    side: at the edge the pool has only just started responding, and a slightly misplaced
    onset would put the wells outside the linear range altogether. Where no window is
    strictly within tolerance, the first admitted by widening that tolerance with the
    standard errors of its slope is used instead, and `strict` records that it was: the
    latitude grows with the noise, so a window admitted this way is weaker evidence that
    the pool responds to dilution than a strict one and is reported as such. A pool with
    no window either way raises.

    With `given`, those wells are read whatever the fit says, and the titer comes from the
    window whose middle dilution is theirs. A dilution lies in several overlapping windows whose
    titers differ, so the window is named rather than left to whichever is found first, and
    naming the middle is the same rule the automatic path applies in reverse. Where no
    strict window centers on their dilution, the titer is not reported.

    Args:
        `windows` (list)
            The windows, from `fit_dilution_windows`.
        `wells` (pandas.DataFrame)
            The QC-passing wells, with ``plate``, ``well`` and ``dilution_factor``
            columns.
        `given` (list or None)
            The configured wells as ``plate:well`` strings, all at one dilution, or `None`
            to locate them. A well is named with its plate because the plates of a pool are
            replicates of one dilution series and so name their wells alike.

    Returns:
        `AnalysisWells`.

    """
    strict = [w for w in windows if w["within_tol"]]
    if given is None:
        # Strict tolerance is preferred over the widening by the error on the slope, as a
        # window admitted only by its error bars is weaker evidence that the pool responds
        # to dilution; it places the range only where nothing else can.
        qualifying = strict or [w for w in windows if w["within_tol_widened"]]
        if not qualifying:
            raise ValueError(
                "no window of dilutions has a slope within tolerance of the ideal "
                f"{IDEAL_SLOPE:+.0f}, even once that tolerance is widened by the standard "
                "errors of the slope, so the pool never clearly responds to being diluted "
                "and there is no linear range to read its composition from. Set "
                "'linear_range_wells' to read it from particular wells anyway"
            )
        window = qualifying[0]["window"]
        dilution = window[len(window) // 2]
        at_dilution = wells[wells["dilution_factor"] == dilution]
        return AnalysisWells(
            wells=sorted(
                at_dilution[["plate", "well"]].itertuples(index=False, name=None)
            ),
            dilution=dilution,
            titer_window=window,
            strict=bool(strict),
            wells_given=False,
        )

    if not given:
        raise ValueError("'linear_range_wells' must name at least one well")
    if len(set(given)) != len(given):
        raise ValueError(
            f"'linear_range_wells' names a well more than once: {sorted(given)}"
        )
    malformed = sorted(w for w in given if len(str(w).split(":")) != 2)
    if malformed:
        raise ValueError(
            f"'linear_range_wells' names {malformed}, but a well is named as "
            "'plate:well', the plates of a pool naming their wells alike"
        )
    given_pairs = [tuple(str(w).split(":")) for w in given]
    known = set(wells[["plate", "well"]].itertuples(index=False, name=None))
    unknown = sorted(set(given_pairs) - known)
    if unknown:
        raise ValueError(
            f"'linear_range_wells' names {[f'{p}:{w}' for (p, w) in unknown]}, which are "
            "not wells of the pool's plates that passed QC; a well that failed QC is not "
            "trusted to measure the pool"
        )
    is_given = pandas.MultiIndex.from_frame(wells[["plate", "well"]]).isin(
        set(given_pairs)
    )
    dilutions = sorted(wells.loc[is_given, "dilution_factor"].unique())
    if len(dilutions) != 1:
        raise ValueError(
            f"the wells of 'linear_range_wells' are at the dilutions {dilutions}, but "
            "the composition is read at one dilution, so they must share one"
        )
    (dilution,) = dilutions
    centered = [
        w["window"] for w in strict if w["window"][len(w["window"]) // 2] == dilution
    ]
    return AnalysisWells(
        wells=sorted(given_pairs),
        dilution=dilution,
        titer_window=centered[0] if centered else None,
        strict=bool(centered),
        wells_given=True,
    )


def pool_samples(samples, pool_wells):
    """The rows of a plate's samples CSV that are this pool's wells.

    One physical plate can hold several pools, each diluted down its own series, and the
    analysis reads the wells it is given as a single series. `pool_wells` picks this pool's
    wells out by whatever column of the samples CSV names the pool.

    Args:
        `samples` (pandas.DataFrame)
            The plate's samples CSV, with a ``well`` column.
        `pool_wells` (dict)
            Maps a samples-CSV column to the value this pool's wells have in it. Empty
            takes the whole plate, which is the case of one pool per plate.

    Returns:
        The rows of `samples` that are this pool's wells, in their original order.

    """
    selected = pandas.Series(True, index=samples.index)
    for column, value in dict(pool_wells).items():
        if column not in samples.columns:
            raise ValueError(
                f"'pool_wells' selects on {column!r}, which is not a column of the "
                f"plate's samples CSV, whose columns are {sorted(samples.columns)}"
            )
        selected &= samples[column].astype(str) == str(value)
    if not selected.any():
        raise ValueError(
            f"'pool_wells' {dict(pool_wells)} matches no well of the plate's samples CSV"
        )
    return samples[selected]


def member_barcodes(viral_library, pool_members):
    """The barcodes of the strains the pool was made from.

    Which strains went in is a fact about how the pool was made and cannot be read off the
    counts: a member that failed to grow and a non-member both come back with no reads, and
    only the first is a low-titer result worth reporting. Selected by a column of the viral
    library rather than listed, so that adding a strain of a subtype does not mean editing
    a list.

    Args:
        `viral_library` (pandas.DataFrame)
            The viral library, with ``barcode`` and ``strain`` columns.
        `pool_members` (dict)
            Maps a library column to the value a member has in it. Empty takes the whole
            library, which is the case of a pool of everything.

    Returns:
        List of the member barcodes.

    A selection that takes only some of a strain's barcodes raises: members are strains,
    and a strain partly in the pool has no meaningful even share.

    """
    selected = pandas.Series(True, index=viral_library.index)
    for column, value in dict(pool_members).items():
        if column not in viral_library.columns:
            raise ValueError(
                f"'pool_members' selects on {column!r}, which is not a column of the "
                f"viral library, whose columns are {sorted(viral_library.columns)}"
            )
        selected &= viral_library[column].astype(str) == str(value)
    if not selected.any():
        raise ValueError("'pool_members' selects no strain of the viral library")
    member_strains = set(viral_library.loc[selected, "strain"])
    partly = sorted(
        strain
        for strain in member_strains
        if not selected[viral_library["strain"] == strain].all()
    )
    if partly:
        raise ValueError(
            f"'pool_members' selects only some of the barcodes of {partly}; a pool is "
            "made from whole strains, and a strain partly in it has no even share"
        )
    return viral_library.loc[selected, "barcode"].tolist()


def pool_units(viral_library, pool_by):
    """How many `pool_by` units the pool holds, which an even share is one over.

    `viral_library` should already be narrowed to the pool's members.

    """
    if pool_by not in POOL_BY_OPTIONS:
        raise ValueError(
            f"'pool_by' must be one of {list(POOL_BY_OPTIONS)}, not {pool_by!r}"
        )
    if pool_by == "stock":
        check_stocks(viral_library)
    return int(viral_library[pool_by].nunique())


def check_stocks(viral_library):
    """Check the ``stock`` column that analyzing or reporting by stock needs.

    A stock is one tube to pipette, so every barcode has to name the tube it came from, and
    a tube has to hold one strain: a stock spanning two strains has no strain share to be
    divided among, and its measured and divided shares could not be compared.

    """
    if "stock" not in viral_library.columns:
        raise ValueError(
            "the viral library needs a 'stock' column giving the tube of every barcode "
            "to analyze or report by stock"
        )
    blank = viral_library["stock"].isna() | (
        viral_library["stock"].astype(str).str.strip() == ""
    )
    if blank.any():
        raise ValueError(
            "the viral library gives no stock for the barcodes "
            f"{sorted(viral_library.loc[blank, 'barcode'])}"
        )
    spanning = sorted(
        stock
        for stock, strains in viral_library.groupby("stock")["strain"]
        if strains.nunique() > 1
    )
    if spanning:
        raise ValueError(
            f"the stocks {spanning} each hold barcodes of more than one strain, but a "
            "stock is one tube and so holds one strain"
        )


def barcode_fractions(counts, members, plate_wells):
    """Each member barcode's mean share of the member reads of the analysis wells.

    The finest unit the pool can be measured in, from which a share per strain, per stock
    or per barcode all follow by summing. Reads of a library strain the pool was not made
    from are contamination rather than composition, so the shares are over member reads
    only: left in, a non-member would dilute every member's share of the pool.

    Args:
        `counts` (pandas.DataFrame)
            Barcode counts, with ``plate``, ``well``, ``barcode`` and ``count`` columns.
        `members` (list)
            The member barcodes, including any with no reads.
        `plate_wells` (list)
            The ``(plate, well)`` of the analysis wells.

    Returns:
        `pandas.DataFrame` with ``barcode`` and ``mean_fraction``, one row per member,
        whose fractions sum to one.

    """
    wells = set(map(tuple, plate_wells))
    in_wells = counts[
        counts["barcode"].isin(members)
        & pandas.MultiIndex.from_frame(counts[["plate", "well"]]).isin(wells)
    ]
    per_well = []
    for key in sorted(wells):
        well_counts = in_wells[
            (in_wells["plate"] == key[0]) & (in_wells["well"] == key[1])
        ]
        total = well_counts["count"].sum()
        if total <= 0:
            raise ValueError(
                f"well {key[1]} of {key[0]} has no reads of any member of the pool, so it "
                "measures nothing about the pool's composition"
            )
        per_well.append(
            well_counts.set_index("barcode")["count"].reindex(members, fill_value=0)
            / total
        )
    return (
        pandas.concat(per_well, axis=1)
        .mean(axis=1)
        .rename("mean_fraction")
        .rename_axis("barcode")
        .reset_index()
    )


def unit_shares(fractions, viral_library, pool_by, equal_volume=None):
    """Each `pool_by` unit's share of the pool, summed over its barcodes.

    Under ``stock`` two shares are reported, because they are two different quantities. The
    measured share is the sum over that stock's own barcodes, which is what the tube
    actually holds. The divided share is the strain's share split equally among its stocks,
    which is what a project balancing per strain pipettes; it holds only where the strain's
    tubes went into the pool at equal volumes, so it is reported only when `equal_volume`
    says they did.

    Args:
        `fractions` (pandas.DataFrame)
            Per-barcode shares, from `barcode_fractions`.
        `viral_library` (pandas.DataFrame)
            The members of the pool, with ``barcode``, ``strain`` and any ``stock``.
        `pool_by` (str)
            The unit to report, one of `POOL_BY_OPTIONS`.
        `equal_volume` (bool or None)
            Whether a strain's stocks went in at equal volumes, which the divided share
            assumes. Required to report one.

    Returns:
        `pandas.DataFrame` with ``unit``, ``strain``, ``measured_share``, and
        ``divided_share`` where it applies.

    """
    if pool_by not in POOL_BY_OPTIONS:
        raise ValueError(
            f"'pool_by' must be one of {list(POOL_BY_OPTIONS)}, not {pool_by!r}"
        )
    columns = ["barcode", "strain"] + (["stock"] if pool_by == "stock" else [])
    labeled = fractions.merge(
        viral_library[columns], on="barcode", validate="one_to_one"
    )
    # The unit is grouped on under the name `unit`, so that a `pool_by` of ``strain`` does
    # not group on the very column the aggregate below also carries: every `pool_by` then
    # gives the same columns, with `strain` always present to look a unit's strain up by.
    shares = (
        labeled.assign(unit=labeled[pool_by])
        .groupby("unit", as_index=False)
        .aggregate(
            measured_share=pandas.NamedAgg("mean_fraction", "sum"),
            strain=pandas.NamedAgg("strain", "first"),
        )
    )
    if pool_by == "stock" and equal_volume:
        per_strain = labeled.groupby("strain")["mean_fraction"].sum()
        stocks_per_strain = labeled.groupby("strain")["stock"].nunique()
        shares["divided_share"] = shares["strain"].map(per_strain / stocks_per_strain)
    elif pool_by != "stock" and equal_volume:
        # a strain or a barcode is its own unit, so there is nothing to divide
        shares["divided_share"] = shares["measured_share"]
    return shares.sort_values("measured_share", ascending=False).reset_index(drop=True)


def diverging_stocks(shares, viral_library, factor):
    """Strains whose stocks' measured shares differ by more than `factor`.

    Two tubes of one strain should hold about the same amount of virus where they went in
    at equal volumes, so a strain whose tubes disagree is worth seeing: it says one tube
    grew poorly, which is what these plates exist to catch. A tube with no reads at all is
    reported without a ratio, there being nothing to divide by; a strain with no reads from
    any of its tubes is a low-titer result rather than a disagreement, and is left out.

    Args:
        `shares` (pandas.DataFrame)
            Per-stock shares, from `unit_shares`, with ``unit``, ``strain`` and
            ``measured_share``.
        `viral_library` (pandas.DataFrame)
            The members of the pool; unused except to require the ``stock`` column, so
            that calling this without one fails here rather than silently reporting
            nothing.
        `factor` (float)
            How many fold the stocks of a strain may differ before being reported.

    Returns:
        `pandas.DataFrame` with ``strain``, ``n_stocks``, ``min_share``, ``max_share`` and
        ``ratio``, which is null where the weakest stock has no reads.

    """
    check_stocks(viral_library)
    rows = []
    for strain, group in shares.groupby("strain"):
        if len(group) < 2 or group["measured_share"].max() <= 0:
            continue
        low = group["measured_share"].min()
        high = group["measured_share"].max()
        ratio = high / low if low > 0 else None
        if ratio is None or ratio > factor:
            rows.append(
                {
                    "strain": strain,
                    "n_stocks": len(group),
                    "min_share": low,
                    "max_share": high,
                    "ratio": ratio,
                }
            )
    return pandas.DataFrame(
        rows, columns=["strain", "n_stocks", "min_share", "max_share", "ratio"]
    ).sort_values("strain", ignore_index=True)
