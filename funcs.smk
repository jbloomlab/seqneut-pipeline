"""Functions used by Snakemake pipeline ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

"""


import copy
import functools
import os


def process_plate(plate, plate_params):
    """Process a plot from the configuration."""

    # Process plate parameters
    req_plate_params = {
        "date",
        "viral_library",
        "neut_standard_set",
        "samples_csv",
        "manual_drops",
        "qc_thresholds",
        "curvefit_params",
        "curvefit_qc",
    }
    if not req_plate_params.issubset(plate_params):
        raise ValueError(f"{plate=} {plate_params=} lacks {req_plate_params=}")
    if plate_params["viral_library"] not in viral_libraries:
        raise ValueError(
            f"{plate=} {plate_params['viral_library']=} not in {viral_libraries=}"
        )
    if plate_params["neut_standard_set"] not in neut_standard_sets:
        raise ValueError(
            f"{plate=} {plate_params['neut_standard_set']=} not in {neut_standard_sets=}"
        )
    plate_d = copy.deepcopy(plate_params)
    plate_d["date"] = str(plate_d["date"])
    if not re.fullmatch("\d{4}\-\d{2}\-\d{2}", str(plate_d["date"])):
        raise ValueError(f"{plate=} {plate_d['date']=} not in YYYY-MM-DD format")

    # Process samples_csv to create the sample data frame
    req_sample_cols = ["well", "serum", "dilution_factor", "replicate", "fastq"]
    samples_df = pd.read_csv(plate_params["samples_csv"], comment="#")
    if not set(req_sample_cols).issubset(samples_df.columns):
        raise ValueError(f"{plate=} {samples_df.columns=} lacks {req_sample_cols=}")

    if samples_df["serum"].isnull().any():
        raise ValueError(f"{plate=} 'samples_csv' has null values in 'serum' column")

    # try to turn columns of ints and NAs into Int64 to avoid ints appearing as flaots
    for col in ["replicate", "dilution_factor"]:
        try:
            samples_df[col] = samples_df[col].astype("Int64")
        except TypeError:
            pass

    # make serum_replicate that defines serum and replicate if needed
    samples_df = (
        samples_df[req_sample_cols]
        .assign(
            one_serum_replicate=lambda x: (
                x.groupby("serum")["replicate"].transform("nunique", dropna=False) == 1
            ),
            serum_replicate=lambda x: x.apply(
                lambda row: (
                    str(row["serum"])
                    + (
                        ""
                        if row["one_serum_replicate"] == 1
                        else f"-{row['replicate']}"
                    )
                ),
                axis=1,
            ),
            sample_noplate=lambda x: x.apply(
                lambda row: (
                    row["serum_replicate"]
                    + (
                        ""
                        if pd.isnull(row["dilution_factor"])
                        else f"_{row['dilution_factor']}"
                    )
                ),
                axis=1,
            ),
            sample=lambda x: x["sample_noplate"].map(lambda s: f"{plate}_{s}"),
            plate=plate,
            plate_replicate=lambda x: x.apply(
                lambda row: (
                    plate
                    + ("" if row["one_serum_replicate"] else f"{-row['replicate']}")
                ),
                axis=1,
            ),
        )
        .drop(columns="one_serum_replicate")
    )

    assert len(samples_df) == samples_df["sample"].nunique(), plate

    # make sure serum_replicate and dilution_factor are unique
    dup_rows = (
        samples_df.assign(
            duplicates=lambda x: (
                x.groupby(["serum_replicate", "dilution_factor"], dropna=False)[
                    "sample"
                ].transform("count")
            ),
        )
        .query("duplicates > 1")
        .drop(columns="duplicates")
    )
    if len(dup_rows):
        raise ValueError(f"{plate=} has duplicated serum / replicates:\n{dup_rows}")

    # make sure dilution_factor is valid
    if not (
        (samples_df["dilution_factor"] >= 1) | (samples_df["serum"] == "none")
    ).all():
        raise ValueError(f"{plate=} has dilution factors not >= 1 for non-none serum")

    # make sure there is at least one "none" sample
    if "none" not in set(samples_df["serum"]):
        raise ValueError(f"{plate=} has no samples with serum set to 'none'")

    # make sure fastqs are unique
    dup_fastqs = (
        samples_df.assign(
            duplicates=lambda x: x.groupby("fastq")["fastq"].transform("count")
        )
        .query("duplicates > 1")
        .drop(columns="duplicates")
    )
    if len(dup_fastqs):
        raise ValueError(f"{plate=} has duplicate FASTQs:\n{dup_fastqs}")

    plate_d["samples"] = samples_df

    return plate_d


@functools.lru_cache
def sera_plates():
    """Get dict keyed by serum with values lists of plates with titers for serum."""
    csv_file = checkpoints.sera_by_plate.get().output.csv
    return (
        pd.read_csv(csv_file)
        .assign(plates=lambda x: x["plates"].str.split(";"))
        .set_index("serum")["plates"]
        .to_dict()
    )
