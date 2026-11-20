"""Make the synthetic pool plates of the test example.

Writes, from one plan of what each pool holds:

  - `data/viral_libraries/pool_library.csv`, the library the pools are made from
  - `fastqs/pool/<plate>/<well>.fastq.gz`, the reads of each well
  - `data/miscellaneous_plates/pool_plate_<n>.csv`, the samples CSV of each plate
  - `expected_pool_composition.csv` and `expected_pool_summary.csv`, what the analyses
    should find

Each plate plants one thing the analysis has to get right. Plate 1 is a pool of the whole
library over a series that saturates at its concentrated end and then responds to dilution,
with reads of a discarded barcode, of misread barcodes, and of barcodes matching nothing, and
one member that did not grow. Plate 2 is a subpool selected by a column, with reads of
non-members, a well too shallow to trust, and a strain whose two stocks differ severalfold.
Plate 3 never clearly responds to dilution, so it is read only because its wells are named.

The library is written here rather than reusing `pdmH1N1_lib2023_loes.csv`, which has no
column `pool_members` could select a whole subpool on and no per-strain annotation for
`addtl_viral_cols_to_show`.

Every read is built rather than sampled, as the reverse complement of its barcode followed by
the exact flanks at full quality, so that every planted read parses and the counts the
pipeline finds are exactly the counts planted here. The output is deterministic, so this can
be rerun to regenerate the files byte for byte. Every path is relative to `test_example`,
which it runs in wherever it is run from:

    python fastqs/pool/make_pool_fastqs.py

"""

import os
import random
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent))

from synthetic_reads import (
    mutate,
    random_barcode,
    read_of,
    split_evenly,
    write_fastq,
)

SEED = 2

LIBRARY_CSV = "data/viral_libraries/pool_library.csv"
NEUT_STANDARD_CSV = "data/neut_standard_sets/loes2023_neut_standards.csv"
SAMPLES_DIR = Path("data/miscellaneous_plates")
FASTQ_DIR = Path("fastqs/pool")
EXPECTED_COMPOSITION_CSV = "expected_pool_composition.csv"
EXPECTED_SUMMARY_CSV = "expected_pool_summary.csv"
CONFIG_YML = "config.yml"

BARCODE_LENGTH = 16

# The library: a strain per entry, with how many stocks it is held as and which subgroup it
# belongs to. Few enough strains that a well clears the default 500 average counts per
# barcode without needing tens of thousands of reads, and two strains of two stocks so that
# the measured and divided per-stock shares differ.
STRAINS = [
    # (strain, n_stocks, subgroup, passage)
    ("A/alpha/1/2020", 1, "group_a", "cell"),
    ("A/alpha/2/2020", 2, "group_a", "cell"),
    ("A/alpha/3/2020", 1, "group_a", "egg"),
    ("A/beta/1/2021", 1, "group_b", "cell"),
    ("A/beta/2/2021", 2, "group_b", "egg"),
    ("A/beta/3/2021", 1, "group_b", "cell"),
]

# Reads of the neutralization standard in a well, spread over its barcodes. Fixed across the
# series, as the standard is the same spike-in in every well; the viral reads are what fall.
NEUT_READS = 3000

# Viral reads of the most concentrated well of a series. A well holds this many until the
# pool starts responding to dilution, after which they fall with the dilution factor, which
# is what puts a known linear range in the data.
#
# It has to leave every well of the linear range averaging clear of the default
# `min_avg_barcode_count_per_well` of 500 over all 18 barcodes counted, the 10 of the
# standard and the 8 of the library, or the range has no run of wells that all pass QC and
# the analysis rightly refuses to read the pool. The margin is worth keeping wide: at 24000
# the well three dilutions past the onset averaged 499.94, which failed by a sixteenth of a
# read and cost the range its last dilution.
VIRAL_READS_AT_PLATEAU = 48000

# the dilutions of a plate's series, and the dilution at which each plate's pool starts
# responding. A pool that never responds has no onset.
DILUTIONS = [2, 4, 8, 16, 32]
REPLICATES = 2


def stock_name(strain, i):
    """The name of stock `i` of `strain`."""
    return f"{strain}_stock{i + 1}"


def build_library(rng):
    """The viral library, one row per barcode, one barcode per stock."""
    rows = []
    for strain, n_stocks, subgroup, passage in STRAINS:
        for i in range(n_stocks):
            rows.append(
                {
                    "barcode": random_barcode(
                        rng,
                        BARCODE_LENGTH,
                        [r["barcode"] for r in rows],
                        min_distance=4,
                    ),
                    "strain": strain,
                    "stock": stock_name(strain, i),
                    "subgroup": subgroup,
                    "passage_history": passage,
                }
            )
    return pd.DataFrame(rows)


def viral_reads_at(dilution, onset):
    """Viral reads of a well at `dilution`, given where the pool starts responding.

    Above `onset` the pool is too concentrated for its counts to fall, so they sit at the
    plateau; from `onset` on they fall with the dilution, which is the linear range. A pool
    with no `onset` never responds and sits at the plateau throughout.

    """
    if onset is None or dilution <= onset:
        return VIRAL_READS_AT_PLATEAU
    return round(VIRAL_READS_AT_PLATEAU * onset / dilution)


def well_name(i_dilution, i_replicate):
    """The well of replicate `i_replicate` at dilution index `i_dilution`."""
    return f"{'ABCDEFGH'[i_replicate]}{i_dilution + 1}"


def plan_plate(plate, library, members, shares, onset, extras):
    """What every well of one plate holds.

    Args:
        `plate` (str)
            The plate's name.
        `library` (pandas.DataFrame)
            The viral library.
        `members` (list)
            The barcodes the pool was made from.
        `shares` (dict)
            Maps a member barcode to its share of the pool's viral reads.
        `onset` (int or None)
            The dilution at which the pool starts responding, or `None` for a pool that
            never does.
        `extras` (dict)
            Per-well extras, mapping a well to a dict of barcode to reads to add on top of
            the pool's own: a non-member, a discarded barcode, a misread, or an unknown
            barcode. A well named here with ``shallow`` instead is given few reads.

    Returns:
        List of dicts, one per well, with its reads by barcode and its dilution.

    """
    neut_barcodes = pd.read_csv(NEUT_STANDARD_CSV)["barcode"].tolist()
    wells = []
    for i_dilution, dilution in enumerate(DILUTIONS):
        for i_replicate in range(REPLICATES):
            well = well_name(i_dilution, i_replicate)
            extra = dict(extras.get(well, {}))
            shallow = extra.pop("shallow", False)
            scale = 0.02 if shallow else 1.0
            viral = round(viral_reads_at(dilution, onset) * scale)
            reads = {}
            # the pool's own barcodes, in the planted proportions
            for barcode, share in shares.items():
                reads[barcode] = round(viral * share)
            # the standard, the same in every well of the series
            for barcode, n in zip(
                neut_barcodes,
                split_evenly(round(NEUT_READS * scale), len(neut_barcodes)),
                strict=True,
            ):
                reads[barcode] = n
            for barcode, n in extra.items():
                reads[barcode] = reads.get(barcode, 0) + n
            wells.append(
                {
                    "plate": plate,
                    "well": well,
                    "dilution_factor": dilution,
                    "replicate": i_replicate + 1,
                    "reads": {b: n for b, n in reads.items() if n > 0},
                }
            )
    return wells


def write_plate(plate, wells, unparsed_barcode, pool_column=None):
    """Write a plate's FASTQs and samples CSV, and return the samples rows.

    `pool_column` names the pool each well holds, for a plate whose analysis selects its
    wells with `pool_wells` rather than taking the whole plate.
    """
    rows = []
    for w in wells:
        sequences = []
        for barcode, n in sorted(w["reads"].items()):
            sequences += [read_of(barcode)] * n
        # a few reads that cannot parse, so the report has more than one fate to draw
        sequences += [read_of(unparsed_barcode, outer="AAAAAA")] * 15
        fastq = FASTQ_DIR / plate / f"{w['well']}.fastq.gz"
        write_fastq(fastq, f"{plate}:{w['well']}", sequences)
        rows.append(
            {
                "well": w["well"],
                "fastq": str(fastq),
                "dilution_factor": w["dilution_factor"],
                "replicate": w["replicate"],
                **({"pool": pool_column} if pool_column else {}),
            }
        )
    pd.DataFrame(rows).to_csv(SAMPLES_DIR / f"{plate}.csv", index=False)
    return rows


def expected_composition(
    plate, wells, library, members, analysis_dilution, pool_by, equal_volume
):
    """What the analysis should find for each member barcode of a pool.

    The per-stock shares are planted only where the analysis reports them, which is when it
    measures the pool by stock: a share per tube says nothing when the unit being balanced
    is the strain or the barcode.

    """
    analysis = [w for w in wells if w["dilution_factor"] == analysis_dilution]
    member_set = set(members)
    fractions = {}
    for barcode in members:
        per_well = []
        for w in analysis:
            member_reads = sum(n for b, n in w["reads"].items() if b in member_set)
            per_well.append(w["reads"].get(barcode, 0) / member_reads)
        fractions[barcode] = sum(per_well) / len(per_well)
    lib = library.set_index("barcode")
    rows = []
    for barcode, fraction in fractions.items():
        stock = lib.loc[barcode, "stock"]
        strain = lib.loc[barcode, "strain"]
        same_stock = [b for b in members if lib.loc[b, "stock"] == stock]
        same_strain = [b for b in members if lib.loc[b, "strain"] == strain]
        n_stocks = len({lib.loc[b, "stock"] for b in same_strain})
        row = {
            "pool": plate,
            "barcode": barcode,
            "strain": strain,
            "stock": stock,
            "mean_fraction": fraction,
        }
        if pool_by == "stock":
            row["stock_measured_share"] = sum(fractions[b] for b in same_stock)
            if equal_volume:
                row["stock_divided_share"] = (
                    sum(fractions[b] for b in same_strain) / n_stocks
                )
        rows.append(row)
    return rows


def main():
    """Write the library, the plates, and the expected answers."""
    os.chdir(Path(__file__).parent.parent.parent)
    rng = random.Random(SEED)

    library = build_library(rng)
    Path(LIBRARY_CSV).parent.mkdir(parents=True, exist_ok=True)
    library.to_csv(LIBRARY_CSV, index=False)

    neut_barcodes = pd.read_csv(NEUT_STANDARD_CSV)["barcode"].tolist()
    known = library["barcode"].tolist() + neut_barcodes
    # the discarded barcode of `config.yml`, read rather than invented so that its planted
    # reads are counted as discarded rather than as matching nothing
    with open(CONFIG_YML) as f:
        discarded = min(yaml.safe_load(f)["discarded_barcodes"])
    unknown = random_barcode(rng, BARCODE_LENGTH, known + [discarded], min_distance=4)
    unparsed = random_barcode(rng, BARCODE_LENGTH, known + [discarded], min_distance=4)

    all_barcodes = library["barcode"].tolist()
    composition_rows = []
    summary_rows = []

    # --- plate 1: the whole library, saturating and then responding ----------------------
    # one member with no reads at all, to check it is reported at zero rather than dropped
    plate_1_members = all_barcodes
    plate_1_shares = {
        b: (0.0 if i == len(all_barcodes) - 1 else 1.0)
        for i, b in enumerate(all_barcodes)
    }
    total = sum(plate_1_shares.values())
    plate_1_shares = {b: s / total for b, s in plate_1_shares.items()}
    onset_1 = 4
    extras_1 = {
        # the discarded spike-in, in every well
        **{
            well_name(i, r): {discarded: 400}
            for i in range(len(DILUTIONS))
            for r in range(REPLICATES)
        },
    }
    # a misread of a member, and a barcode matching nothing, in the wells read
    analysis_dilution_1 = DILUTIONS[DILUTIONS.index(onset_1) + 1]
    for r in range(REPLICATES):
        well = well_name(DILUTIONS.index(analysis_dilution_1), r)
        extras_1[well] = {
            **extras_1[well],
            mutate(all_barcodes[0], 3): 120,
            unknown: 200,
        }
    wells_1 = plan_plate(
        "pool_plate_1", library, plate_1_members, plate_1_shares, onset_1, extras_1
    )
    write_plate("pool_plate_1", wells_1, unparsed)
    composition_rows += expected_composition(
        "pool_plate_1",
        wells_1,
        library,
        plate_1_members,
        analysis_dilution_1,
        "strain",
        True,
    )
    summary_rows.append(
        {
            "pool": "pool_plate_1",
            "analysis_dilution": analysis_dilution_1,
            "linear_range_dilutions": " ".join(
                str(d)
                for d in DILUTIONS[
                    DILUTIONS.index(analysis_dilution_1)
                    - 1 : DILUTIONS.index(analysis_dilution_1)
                    + 2
                ]
            ),
            "linear_range_strict": True,
            "analysis_wells_given": False,
            "n_strains": library["strain"].nunique(),
        }
    )

    # --- plate 2: a subpool, with non-members, a shallow well, and uneven stocks ----------
    members_2 = library.loc[library["subgroup"] == "group_a", "barcode"].tolist()
    non_members = [b for b in all_barcodes if b not in members_2]
    # the two stocks of the two-stock strain differ severalfold, for the divergence warning
    two_stock_barcodes = library.loc[
        (library["subgroup"] == "group_a")
        & library["strain"].isin(
            library.groupby("strain")["stock"].nunique().pipe(lambda s: s[s > 1]).index
        ),
        "barcode",
    ].tolist()
    plate_2_shares = {}
    for b in members_2:
        if b == two_stock_barcodes[0]:
            plate_2_shares[b] = 4.0
        elif b in two_stock_barcodes:
            plate_2_shares[b] = 0.5
        else:
            plate_2_shares[b] = 1.0
    total = sum(plate_2_shares.values())
    plate_2_shares = {b: s / total for b, s in plate_2_shares.items()}
    onset_2 = 4
    analysis_dilution_2 = DILUTIONS[DILUTIONS.index(onset_2) + 1]
    extras_2 = {
        well_name(i, r): {non_members[0]: 300}
        for i in range(len(DILUTIONS))
        for r in range(REPLICATES)
    }
    # one replicate of the most dilute point is too shallow to trust
    extras_2[well_name(len(DILUTIONS) - 1, 1)] = {
        **extras_2[well_name(len(DILUTIONS) - 1, 1)],
        "shallow": True,
    }
    wells_2 = plan_plate(
        "pool_plate_2", library, members_2, plate_2_shares, onset_2, extras_2
    )
    # named in a `pool` column so that this analysis selects its wells with `pool_wells`,
    # as a plate holding several pools must
    write_plate("pool_plate_2", wells_2, unparsed, pool_column="subgroup_pool")
    composition_rows += expected_composition(
        "pool_plate_2", wells_2, library, members_2, analysis_dilution_2, "stock", True
    )
    summary_rows.append(
        {
            "pool": "pool_plate_2",
            "analysis_dilution": analysis_dilution_2,
            "linear_range_dilutions": " ".join(
                str(d)
                for d in DILUTIONS[
                    DILUTIONS.index(analysis_dilution_2)
                    - 1 : DILUTIONS.index(analysis_dilution_2)
                    + 2
                ]
            ),
            "linear_range_strict": True,
            "analysis_wells_given": False,
            "n_strains": library.loc[
                library["subgroup"] == "group_a", "strain"
            ].nunique(),
        }
    )

    # --- plate 3: never responds to dilution, so its wells must be named -----------------
    plate_3_shares = {b: 1 / len(all_barcodes) for b in all_barcodes}
    wells_3 = plan_plate(
        "pool_plate_3", library, all_barcodes, plate_3_shares, None, {}
    )
    write_plate("pool_plate_3", wells_3, unparsed)
    given_dilution = DILUTIONS[2]
    composition_rows += expected_composition(
        "pool_plate_3", wells_3, library, all_barcodes, given_dilution, "barcode", False
    )
    summary_rows.append(
        {
            "pool": "pool_plate_3",
            "analysis_dilution": given_dilution,
            # no run of a flat series is within tolerance, so there is no titer to report
            "linear_range_dilutions": None,
            "linear_range_strict": False,
            "analysis_wells_given": True,
            "n_strains": library["strain"].nunique(),
        }
    )

    pd.DataFrame(composition_rows).to_csv(
        EXPECTED_COMPOSITION_CSV, index=False, float_format="%.4g"
    )
    pd.DataFrame(summary_rows).to_csv(EXPECTED_SUMMARY_CSV, index=False)

    print(f"wrote {LIBRARY_CSV} with {len(library)} barcodes")
    print(
        f"the wells of plate 3 to name in `linear_range_wells`: "
        f"{[well_name(DILUTIONS.index(given_dilution), r) for r in range(REPLICATES)]}"
    )
    for plate, wells in [
        ("pool_plate_1", wells_1),
        ("pool_plate_2", wells_2),
        ("pool_plate_3", wells_3),
    ]:
        reads = sum(sum(w["reads"].values()) for w in wells)
        print(f"wrote {plate}: {len(wells)} wells, {reads} planted reads")


if __name__ == "__main__":
    main()
