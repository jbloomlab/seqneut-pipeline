"""Make the synthetic single-well plates of the test example.

Writes, from one plan of what each well holds:

  - `fastqs/single_well/<plate>/<well>.fastq.gz`, the reads of each well
  - `data/miscellaneous_plates/<plate>.csv`, the samples CSV of each plate
  - `data/viral_libraries/pdmH1N1_lib2023_loes_final.csv`, the library carried forward
  - `expected_single_well_composition.csv`, what the analysis should find in each well

Each well plants one situation the analysis has to get right: spill from a neighboring
well, material from far across the plate, from another plate of the group, or of a strain
in no well, unmatched barcodes near and far from what belongs, a discarded spike-in barcode
and its misreads, an expected barcode with no reads, an empty well, and a shallow one. Every
read is built rather than sampled, as the reverse complement of its barcode followed by the
exact flanks at full quality, so that every planted read parses and the counts the pipeline
finds are exactly the counts planted here. A few reads that cannot parse are added too, so
that the report has more than one fate to draw.

The output is deterministic, so this can be rerun to regenerate the files byte for byte.
Every path, both those it writes and the FASTQ paths in the samples CSVs, is relative to
`test_example`, which it runs in wherever it is run from:

    python fastqs/single_well/make_single_well_fastqs.py

"""

import gzip
import os
import random
from pathlib import Path

import pandas as pd
import yaml

SEED = 1

LIBRARY_CSV = "data/viral_libraries/pdmH1N1_lib2023_loes.csv"
FINAL_LIBRARY_CSV = "data/viral_libraries/pdmH1N1_lib2023_loes_final.csv"
NEUT_STANDARD_CSV = "data/neut_standard_sets/loes2023_neut_standards.csv"
SAMPLES_DIR = Path("data/miscellaneous_plates")
FASTQ_DIR = Path("fastqs/single_well")
EXPECTED_CSV = "expected_single_well_composition.csv"

# The flanks of `illumina_barcode_parser_params` in `config.yml`, with the `upstream2` the
# single-well plates are given there, as they appear in a read: the read is the reverse
# complement of barcode, upstream, upstream2, as `bc_orientation` is R2.
UPSTREAM = "CTCCCTACAATGTCGGATTTGTATTTAATAG"
UPSTREAM2 = "GCTACA"

# reads of the neutralization standard in an ordinary well, spread over its barcodes
NEUT_READS = 400

# reads of a well's own unit in an ordinary well, spread over its expected barcodes
OWN_READS = 3000

# The configuration whose `analyze_single_well_infections` entries give the thresholds each
# plate's wells are flagged by, read here only to write down which wells should be flagged.
CONFIG_YML = "config.yml"
FLAG_THRESHOLDS = [
    "max_other_library_frac",
    "max_unmatched_frac",
    "min_parsed_reads_per_well",
    "min_neut_standard_reads_for_titer",
]


def flag_thresholds_by_plate():
    """The thresholds each plate's wells are flagged by, from its analysis in `CONFIG_YML`."""
    with open(CONFIG_YML) as f:
        analyses = yaml.safe_load(f)["analyze_single_well_infections"]
    return {
        plate: {key: entry[key] for key in FLAG_THRESHOLDS}
        for entry in analyses.values()
        for plate in entry["miscellaneous_plates"]
    }


def revcomp(seq):
    """Reverse complement of `seq`."""
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def hamming(s1, s2):
    """Hamming distance between equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=True))


def mutate(barcode, position):
    """`barcode` with the nucleotide at `position` changed."""
    new = "C" if barcode[position] != "C" else "G"
    return barcode[:position] + new + barcode[position + 1 :]


def random_barcode(rng, length, far_from, min_distance, prefix=""):
    """A random barcode starting `prefix`, at least `min_distance` from all of `far_from`."""
    while True:
        barcode = prefix + "".join(
            rng.choice("ACGT") for _ in range(length - len(prefix))
        )
        if all(hamming(barcode, other) >= min_distance for other in far_from):
            return barcode


def read_of(barcode, outer=UPSTREAM2):
    """The sequence of a read of `barcode`."""
    return revcomp(barcode) + revcomp(UPSTREAM) + revcomp(outer)


def split_evenly(total, n):
    """Split `total` reads as evenly as possible into `n` parts."""
    return [total // n + (1 if i < total % n else 0) for i in range(n)]


def final_library_status(own_barcodes, final_barcodes, applies):
    """How much of a well's unit the final library keeps, if one applies to its plate."""
    if not applies:
        return {"n_expected_barcodes_final": None, "kept_in_final_library": None}
    n_final = len(set(own_barcodes) & final_barcodes)
    if n_final == 0:
        kept = "none"
    elif n_final < len(own_barcodes):
        kept = "some barcodes"
    else:
        kept = "all barcodes"
    return {"n_expected_barcodes_final": n_final, "kept_in_final_library": kept}


def main():
    """Write the plates."""
    rng = random.Random(SEED)
    library = pd.read_csv(LIBRARY_CSV)
    neut_barcodes = pd.read_csv(NEUT_STANDARD_CSV)["barcode"].tolist()
    thresholds_by_plate = flag_thresholds_by_plate()
    strain_barcodes = library.groupby("strain")["barcode"].apply(sorted).to_dict()
    # strains in a fixed order, the ones with three barcodes first so that the wells that
    # need several barcodes have them
    strains = sorted(strain_barcodes, key=lambda s: (-len(strain_barcodes[s]), s))
    s = dict(enumerate(strains, start=1))  # s[1], s[2], ... name the strains used

    whitelist = [*library["barcode"], *neut_barcodes]
    bclen = len(whitelist[0])
    # More than 2 nt from every barcode counted, so that its misread one nt off it is more
    # than 1 nt from every barcode counted, and so an unmatched read far from what belongs:
    # only the exact barcode is discarded. It starts
    # with GG, as the spike-in barcodes that are discarded in practice do: they are the
    # ones that sequence badly beside indices that also hold GGs.
    discarded = random_barcode(rng, bclen, whitelist, 3, prefix="GG")
    discarded_misread = mutate(discarded, 5)
    # an unmatched barcode too far from everything to be put down to sequencing error
    foreign = random_barcode(rng, bclen, [*whitelist, discarded], 3)

    # What each well holds: its unit (None for an empty well) and its reads, as
    # barcode -> count, plus the expected answer for the source of its top contaminant.
    # `own` and `neut` default to the ordinary amounts; `own_split` overrides how the own
    # reads fall over the unit's barcodes.
    plan = {
        "single_well_plate_1": {
            "B2": {"unit": s[1]},  # clean
            "B3": {"unit": s[2], "other": {s[1]: 60}, "origin": "1 well away"},
            "B4": {"unit": s[3], "other": {s[6]: 40}, "origin": "over 2 wells away"},
            "B5": {
                "unit": s[4],
                "near_own": 30,  # own barcodes misread by one nt
                "foreign": 200,  # far from everything: flags the well as unmatched
            },
            "B6": {"unit": s[5], "discarded": 500, "discarded_misread": 50},
            "B7": {
                "unit": None,  # empty: only the standard, and a little of B5's strain
                "other": {s[4]: 20},
                "origin": "2 wells away",
            },
            "D7": {"unit": s[6]},  # home of B4's far spill
            "C2": {"unit": s[7], "other": {s[11]: 50}, "origin": "another plate here"},
            "C3": {"unit": s[8], "other": {s[20]: 30}, "origin": "not assayed here"},
            "C4": {"unit": s[9], "own_split": "first_barcode_absent"},
            "C5": {"unit": s[10], "own": 450, "neut": 50},  # shallow; imprecise titer
        },
        "single_well_plate_2": {
            "B2": {"unit": s[11]},
            "B3": {"unit": s[12], "other": {s[2]: 40}, "origin": "another plate here"},
            # the same strain as plate 1's B2, so its axis label has to name the plate
            "B5": {"unit": s[1]},
        },
    }
    # One barcode per well: two wells of one strain's barcodes, the second holding some of
    # the first's, which is foreign as the well was meant to hold its one barcode only.
    # These two barcodes share a tube, which analyzing by barcode ignores: the well was
    # infected from one barcode and so the other does not belong in it either.
    sibling_a, sibling_b = strain_barcodes[s[13]][:2]
    assert (
        library.loc[library["barcode"] == sibling_a, "stock_id"].item()
        == library.loc[library["barcode"] == sibling_b, "stock_id"].item()
    ), "the barcode plate's two wells of one strain hold barcodes of one tube"
    plan_barcodes = {
        "single_well_plate_3": {
            "B2": {"unit": sibling_a, "strain": s[13]},
            "B3": {
                "unit": sibling_b,
                "strain": s[13],
                "other_barcodes": {sibling_a: 60},
                "origin": "1 well away",
            },
            "B4": {"unit": strain_barcodes[s[14]][0], "strain": s[14]},
        },
    }

    # One tube per well. The library's `stock_id` records which barcodes share a tube: the
    # first two barcodes of s[15] in library order share its first and its third is a
    # second tube, every other barcode being a tube of its own. The well infected from the
    # shared tube holds both of its barcodes, and the wells of s[15] are labeled by which
    # of its tubes they hold. The shared tube is in a second well of the plate too,
    # which is first in the order the plate was filled but farther from the well its
    # spill is found in, so that spill has to be traced to the nearer of the two.
    stock_barcodes = library.groupby("stock_id")["barcode"].apply(list).to_dict()
    in_library_order = library.loc[library["strain"] == s[15], "barcode"].tolist()
    shared_tube, second_tube = (
        library.loc[library["barcode"] == barcode, "stock_id"].item()
        for barcode in [in_library_order[0], in_library_order[2]]
    )
    assert (
        len(stock_barcodes[shared_tube]) == 2
    ), "s[15]'s first tube holds two barcodes"
    plan_stocks = {
        "single_well_plate_4": {
            "B2": {"unit": shared_tube, "strain": s[15]},
            "B3": {
                "unit": second_tube,
                "strain": s[15],
                "other_barcodes": {in_library_order[0]: 60},
                "origin": "1 well away",
            },
            # last, so that adding it left the other wells unchanged
            "A1": {"unit": shared_tube, "strain": s[15]},
        },
    }

    # the library carried forward: s[9] dropped entirely, s[7] with one barcode dropped
    final = library[
        (library["strain"] != s[9]) & (library["barcode"] != strain_barcodes[s[7]][0])
    ]
    final.to_csv(FINAL_LIBRARY_CSV, index=False)
    final_barcodes = set(final["barcode"])
    # the plates analyzed against the final library, as `config.yml` configures them
    final_library_plates = set(plan)

    expected = []
    # the plate recording stocks comes last, so that adding it left the others unchanged
    for plate, wells in [*plan.items(), *plan_barcodes.items(), *plan_stocks.items()]:
        # the column of the samples CSV naming what each well was infected from, which is
        # the viral library's column for the `analyze_by` of the plate's analysis in
        # `CONFIG_YML`
        unit_column = (
            "barcode"
            if plate in plan_barcodes
            else "stock_id" if plate in plan_stocks else "strain"
        )
        thresholds = thresholds_by_plate[plate]
        sample_rows = []
        (FASTQ_DIR / plate).mkdir(parents=True, exist_ok=True)
        for well, w in wells.items():
            unit = w["unit"]
            reads = {}  # barcode -> count

            def add(barcode, n, reads=reads):
                reads[barcode] = reads.get(barcode, 0) + n

            for bc, n in zip(
                neut_barcodes,
                split_evenly(w.get("neut", NEUT_READS), len(neut_barcodes)),
                strict=True,
            ):
                add(bc, n)
            if unit is not None:
                if plate in plan_stocks:
                    own_barcodes = sorted(stock_barcodes[unit])
                elif plate in plan_barcodes:
                    own_barcodes = [unit]
                else:
                    own_barcodes = strain_barcodes[unit]
                own = w.get("own", OWN_READS)
                if w.get("own_split") == "first_barcode_absent":
                    split = [0, *split_evenly(own, len(own_barcodes) - 1)]
                else:
                    split = split_evenly(own, len(own_barcodes))
                for bc, n in zip(own_barcodes, split, strict=True):
                    add(bc, n)
                absent = [
                    bc for bc, n in zip(own_barcodes, split, strict=True) if not n
                ]
            else:
                own_barcodes = []
                absent = []
            for other_strain, n in w.get("other", {}).items():
                for bc, k in zip(
                    strain_barcodes[other_strain],
                    split_evenly(n, len(strain_barcodes[other_strain])),
                    strict=True,
                ):
                    add(bc, k)
            for bc, n in w.get("other_barcodes", {}).items():
                add(bc, n)
            near_own = w.get("near_own", 0)
            if near_own:
                add(mutate(own_barcodes[0], 3), near_own)
            add(foreign, w.get("foreign", 0))
            add(discarded, w.get("discarded", 0))
            add(discarded_misread, w.get("discarded_misread", 0))
            reads = {bc: n for (bc, n) in reads.items() if n > 0}

            sequences = [read_of(bc) for (bc, n) in reads.items() for _ in range(n)]
            # reads that cannot be parsed: a wrong outer flank, and a low-quality barcode
            n_unparsed = 15
            sequences += [read_of(neut_barcodes[0], outer="AAAAAA")] * n_unparsed
            rng.shuffle(sequences)
            qualities = ["F" * len(seq) for seq in sequences]
            sequences.append(read_of(neut_barcodes[0]))
            qualities.append("#" * bclen + "F" * (len(sequences[-1]) - bclen))
            fastq = FASTQ_DIR / plate / f"{well}.fastq.gz"
            # `mtime=0` so that the gzip header, and so the file, is the same every run
            with (
                open(fastq, "wb") as raw,
                gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as f,
            ):
                for i, (seq, qual) in enumerate(zip(sequences, qualities, strict=True)):
                    f.write(f"@{plate}:{well}:{i}\n{seq}\n+\n{qual}\n".encode())

            if unit_column == "strain":
                row = {"well": well, "strain": unit}
            else:
                row = {"well": well, "strain": w["strain"], unit_column: unit}
            sample_rows.append({**row, "fastq": str(fastq), "dilution_factor": 1})

            # the answers, by category
            library_reads = {
                bc: n for bc, n in reads.items() if bc in set(library["barcode"])
            }
            expected_reads = sum(
                n for bc, n in library_reads.items() if bc in own_barcodes
            )
            other_library_reads = sum(library_reads.values()) - expected_reads
            neut_reads = sum(n for bc, n in reads.items() if bc in neut_barcodes)
            discarded_reads = w.get("discarded", 0)
            far = w.get("foreign", 0) + w.get("discarded_misread", 0)
            non_neut = expected_reads + other_library_reads + near_own + far
            parsed = non_neut + neut_reads
            frac_other = other_library_reads / non_neut if non_neut else float("nan")
            frac_unmatched = (near_own + far) / non_neut if non_neut else float("nan")
            expected.append(
                {
                    "plate": plate,
                    "well": well,
                    "n_expected_barcodes": len(own_barcodes),
                    "expected_reads": expected_reads,
                    "neut_standard_reads": neut_reads,
                    "other_library_reads": other_library_reads,
                    "unmatched_near_expected_reads": near_own,
                    "unmatched_far_expected_reads": far,
                    "discarded_reads": discarded_reads,
                    "flagged_other_library": frac_other
                    > thresholds["max_other_library_frac"],
                    "flagged_unmatched": frac_unmatched
                    > thresholds["max_unmatched_frac"],
                    "too_shallow": parsed < thresholds["min_parsed_reads_per_well"],
                    "titer_imprecise": unit is not None
                    and neut_reads < thresholds["min_neut_standard_reads_for_titer"],
                    "titer": expected_reads / neut_reads if unit is not None else None,
                    "top_other_unit_origin": w.get("origin"),
                    # barcodes the well should hold that were given no reads at all
                    "absent_barcodes": ";".join(absent) or None,
                    **final_library_status(
                        own_barcodes,
                        final_barcodes,
                        unit is not None and plate in final_library_plates,
                    ),
                }
            )
        pd.DataFrame(sample_rows).to_csv(SAMPLES_DIR / f"{plate}.csv", index=False)

    pd.DataFrame(expected).to_csv(EXPECTED_CSV, index=False)
    print(f"discarded barcode: {discarded}")


if __name__ == "__main__":
    # `test_example`, two levels up from this script
    os.chdir(Path(__file__).resolve().parents[2])
    main()
