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

# the unit a pool can be balanced by, and so the unit an equal share is an equal share of
POOL_BY_OPTIONS = ("strain", "barcode", "stock")


@dataclass(frozen=True)
class AnalysisWells:
    """The wells a pool's composition is read from, and where its titer comes from.

    Attributes:
        `wells` (list)
            The ``(plate, well)`` of each well to read, all at one dilution. Empty where no
            window of the series is within tolerance, there being no dilution the pool can
            be read at.
        `gap_blocked` (list)
            The windows that close over a dilution the plate holds, and so cannot be
            selected whatever their slope. Reported, as one of them is what the fit would
            have considered had those wells passed QC.
        `dilution` (int or float or None)
            The dilution those wells are at, or `None` where there are none.
        `titer_window` (list or None)
            The dilutions of the window the titer is computed over, or `None` where no
            strict window gives one and the titer is therefore not reported.
        `strict` (bool)
            Whether a window within tolerance gives the titer. True where the wells were
            located from the fit, which only such a window can do, and false where no
            window is within tolerance or where configured wells sit at a dilution none
            centers.
        `wells_given` (bool)
            Whether the wells were configured rather than located from the fit.

    """

    wells: list
    dilution: float
    titer_window: list
    strict: bool
    wells_given: bool
    gap_blocked: tuple = ()


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


def dilution_means(wells):
    """The geometric mean of the neutralization-standard odds at each dilution.

    Replicate wells of a dilution are averaged before the slope and the titer are taken
    from them, so that each dilution carries the same weight however many replicates of it
    were measured: a dilution left out of some replicates, or whose wells mostly failed QC,
    would otherwise pull the fit towards itself.

    The mean is taken of ``log10(odds)``, which makes it the geometric mean of the odds.
    That is the right centre on a log-log fit, being the arithmetic mean on the axis the
    fit is on, and it is also far less moved by a single wild well than the arithmetic mean
    of the odds themselves.

    Args:
        `wells` (pandas.DataFrame)
            Wells with ``dilution_factor`` and ``neut_standard_odds`` columns.

    Returns:
        Dict mapping each dilution to its mean odds, in increasing order of dilution.

    """
    nonpositive = sorted(v for v in wells["neut_standard_odds"] if not v > 0)
    if nonpositive:
        raise ValueError(
            f"the neutralization standard odds {nonpositive} are not positive, and a "
            "logarithm is defined only for a positive value"
        )
    # `tolist` rather than leaving these as `numpy` numbers: a `numpy` dilution carries its
    # type into the reported window and writes as `np.int64(128)`
    means = (
        wells.assign(_log_odds=[math.log10(v) for v in wells["neut_standard_odds"]])
        .groupby("dilution_factor")["_log_odds"]
        .mean()
    )
    return {
        dilution: 10**log_odds
        for dilution, log_odds in zip(means.index.tolist(), means.tolist(), strict=True)
    }


def spans_gap(window, plate_dilutions):
    """Whether `window` skips a dilution the plate holds.

    A window is a run of consecutive dilutions among those that reach the fit, so a
    dilution whose every well failed QC closes the window over the gap. Such a window
    covers a wider span of the series than its length says, and its middle is not the
    middle of the dilutions actually measured. Where the gap sits next to the saturated
    end this is not merely imprecise: straddling it pulls the slope into tolerance, and the
    window being the more concentrated one it is then selected first, so the composition
    would be read from a saturated well, which is the one thing locating the linear range
    exists to prevent.

    Args:
        `window` (list)
            The dilutions of the window.
        `plate_dilutions` (iterable)
            Every dilution the plate holds, whether or not its wells reached the fit.

    Returns:
        True where the plate holds a dilution inside the window's span that the window
        does not cover.

    """
    held = [d for d in plate_dilutions if window[0] <= d <= window[-1]]
    return len(held) > len(window)


def fit_dilution_windows(
    wells,
    window_dilutions=DEFAULT_LINEARITY_WINDOW_DILUTIONS,
    slope_tol=DEFAULT_LINEARITY_SLOPE_TOL,
    plate_dilutions=None,
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
        `plate_dilutions` (iterable or None)
            Every dilution the plate holds, including any whose wells all failed QC. A
            window that closes over one of those is reported but is not within tolerance,
            so the composition is never read from it; see `spans_gap`. `None` takes the
            dilutions of `wells`, under which no window spans a gap.

    Returns:
        List of dicts, one per window from the most concentrated end, each with the window's
        ``window`` of dilutions, its ``slope``, ``slope_se`` and ``intercept``, which
        together place its fitted line, its ``deviation`` from `IDEAL_SLOPE`, whether it
        ``spans_gap``, and whether it is ``within_tol``. Empty where the series is shorter
        than one window.

    """
    if window_dilutions < 3 or window_dilutions % 2 == 0:
        raise ValueError(
            f"'linearity_window_dilutions' must be an odd number of at least 3, not "
            f"{window_dilutions}: fewer than three dilutions cannot give a slope and its "
            "error, and an even window has no middle dilution to read the composition from"
        )
    # Each dilution once, weighted the same however many replicates of it were measured,
    # and in increasing order, so that a window is a run of consecutive dilutions.
    mean_odds = dilution_means(wells)
    dilutions = sorted(mean_odds)
    if plate_dilutions is None:
        plate_dilutions = dilutions
    plate_dilutions = sorted(set(plate_dilutions))
    windows = []
    for i in range(len(dilutions) - window_dilutions + 1):
        window = dilutions[i : i + window_dilutions]
        slope, slope_se, intercept = log_log_slope(
            window, [mean_odds[d] for d in window]
        )
        deviation = abs(slope - IDEAL_SLOPE)
        gapped = spans_gap(window, plate_dilutions)
        windows.append(
            {
                "window": window,
                "slope": slope,
                "slope_se": slope_se,
                # with the slope, where this window's own fitted line is drawn
                "intercept": intercept,
                "deviation": deviation,
                "spans_gap": gapped,
                # A window that closes over a dilution the plate holds is not trusted
                # however close its slope comes: it does not describe the run of dilutions
                # it appears to, and taking it could read the composition from inside
                # saturation.
                "within_tol": deviation <= slope_tol and not gapped,
            }
        )
    return windows


def analysis_wells(windows, wells, given=None):
    """Choose the wells the composition is read from, and the window the titer comes from.

    Without `given`, the wells are the middle dilution of the first strictly linear window,
    rather than its most concentrated edge, so that there is a dilution of margin on either
    side: at the edge the pool has only just started responding, and a slightly misplaced
    onset would put the wells outside the linear range altogether. Where no window is
    within tolerance the pool never clearly responds to being diluted, and no wells are
    returned: there is no dilution at which what it holds can be read, and reading it at
    one anyway would report a composition the series does not support. The caller reports
    that rather than a composition, and `linear_range_wells` names the wells to read it
    from where that is the intent.

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
            to locate them. A well of a pool measured from a single plate may be named on
            its own, as the report names it, there being no second plate to tell it from;
            where the pool has several plates the plate is required, they being replicates
            of one dilution series and so naming their wells alike.

    Returns:
        `AnalysisWells`.

    """
    strict = [w for w in windows if w["within_tol"]]
    # Close to the ideal slope but for closing over a dilution the plate holds, which is
    # worth reporting: such a window is what the fit would have taken had those wells
    # passed QC, and the reader is the one who knows whether the gap matters.
    gap_blocked = tuple(w["window"] for w in windows if w["spans_gap"])
    if given is None:
        if not strict:
            return AnalysisWells(
                wells=[],
                dilution=None,
                titer_window=None,
                strict=False,
                wells_given=False,
                gap_blocked=gap_blocked,
            )
        window = strict[0]["window"]
        dilution = window[len(window) // 2]
        at_dilution = wells[wells["dilution_factor"] == dilution]
        return AnalysisWells(
            wells=sorted(
                at_dilution[["plate", "well"]].itertuples(index=False, name=None)
            ),
            dilution=dilution,
            titer_window=window,
            # located from the fit, which only a window within tolerance can do
            strict=True,
            wells_given=False,
            gap_blocked=gap_blocked,
        )

    if not given:
        raise ValueError("'linear_range_wells' must name at least one well")
    if len(set(given)) != len(given):
        raise ValueError(
            f"'linear_range_wells' names a well more than once: {sorted(given)}"
        )
    # A single-plate pool's wells are reported without a plate, and copying the reported
    # wells back into `linear_range_wells` to pin them is the obvious way to use this, so a
    # bare well name is taken to mean that one plate.
    plates = sorted(set(wells["plate"]))
    one_plate = len(plates) == 1
    malformed = sorted(
        w
        for w in given
        if len(str(w).split(":")) != 2 and not (one_plate and ":" not in str(w))
    )
    if malformed:
        named_as = "'well' or 'plate:well'" if one_plate else "'plate:well'"
        raise ValueError(
            f"'linear_range_wells' names {malformed}, but a well is named as {named_as}, "
            "the plates of a pool naming their wells alike"
        )
    given_pairs = [
        tuple(str(w).split(":")) if ":" in str(w) else (plates[0], str(w))
        for w in given
    ]
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
    # `tolist` rather than `unique`, so the dilutions are Python numbers: a `numpy` one
    # writes into the message as `np.int64(32)`
    dilutions = sorted(set(wells.loc[is_given, "dilution_factor"].tolist()))
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
        gap_blocked=gap_blocked,
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


def unit_shares(fractions, viral_library, pool_by):
    """Each `pool_by` unit's share of the pool, summed over its barcodes.

    Under ``stock`` the share is what that one tube holds, measured from its own barcodes,
    which is the quantity a volume per tube is read off. A strain held as several tubes
    therefore gets a share for each rather than one share to divide between them.

    Args:
        `fractions` (pandas.DataFrame)
            Per-barcode shares, from `barcode_fractions`.
        `viral_library` (pandas.DataFrame)
            The members of the pool, with ``barcode``, ``strain`` and any ``stock``.
        `pool_by` (str)
            The unit to report, one of `POOL_BY_OPTIONS`.

    Returns:
        `pandas.DataFrame` with ``unit``, ``strain`` and ``measured_share``.

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
    return shares.sort_values("measured_share", ascending=False).reset_index(drop=True)


def diverging_stocks(shares, viral_library, factor):
    """Strains whose stocks' measured shares differ by more than `factor`.

    A strain whose tubes disagree is worth seeing whatever volumes they went in at. Added
    at equal volumes they should hold about the same amount of virus, so a disagreement
    says one tube grew poorly, which is what these plates exist to catch; added to a
    chosen ratio they should differ by about that ratio, so a disagreement far wider than
    any ratio intended says the same thing. Which it is cannot be read off the shares, so
    the ratio is reported and left to a reader who knows what went in.

    A tube with no reads at all is reported without a ratio, there being nothing to divide
    by; a strain with no reads from any of its tubes is a low-titer result rather than a
    disagreement, and is left out.

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
        # `>` rather than `>=`: stocks exactly `factor`-fold apart are within what the
        # threshold allows, as `max_stock_share_ratio` names the largest ratio accepted
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
