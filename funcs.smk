"""Functions used by Snakemake pipeline ``seqneut-pipeline`.

Designed to be included in another ``Snakefile`` that specifies the config.

Note: This file relies on imports from seqneut-pipeline.smk (pandas as pd, re).

"""

import copy
import functools
import os


def stringify_plate_dates(config):
    """Convert the ``date`` of each plate in ``config`` to a string.

    YAML parses unquoted dates (e.g. ``2023-08-01``) into ``datetime.date`` objects,
    which ``snakemake`` cannot JSON serialize when it hashes the config at startup.

    """
    for config_key in ["plates", "miscellaneous_plates"]:
        # `or {}` rather than a `get` default, as either key can be present but null
        for plate_d in (config.get(config_key) or {}).values():
            if "date" in plate_d:
                plate_d["date"] = str(plate_d["date"])


def check_no_whitespace(values, description):
    """Raise a ``ValueError`` if any of ``values`` contains whitespace.

    Args:
        `values` (iterable)
            The values to check.
        `description` (str)
            Names these values in the error message.

    """
    with_whitespace = sorted(
        {str(val) for val in values if any(char.isspace() for char in str(val))}
    )
    if with_whitespace:
        raise ValueError(
            f"{description} cannot contain whitespace, as it would break the output "
            f"paths or the QC-drop YAML keys built from them: {with_whitespace}"
        )


def get_collapse_strain_barcodes(config):
    """Get the boolean value of the optional ``collapse_strain_barcodes`` config key.

    Defaults to ``False`` if the key is absent or null. Strings are also accepted,
    because ``snakemake`` only converts the value of ``--config`` to a bool if it is
    capitalized as in Python, so ``--config collapse_strain_barcodes=true`` would
    otherwise be the (always truthy) string ``"true"``.

    """
    val = config.get("collapse_strain_barcodes", False)
    if val is None:
        return False
    if isinstance(val, bool):
        return val
    if val in {"true", "True"}:
        return True
    if val in {"false", "False"}:
        return False
    raise ValueError(f"invalid {val=} for 'collapse_strain_barcodes', must be a bool")


def process_miscellaneous_plates(misc_plates_d):
    """Process the dictionary of miscellaneous_plates."""
    misc_plates = {}
    req_keys = {"viral_library", "neut_standard_set", "samples_csv"}
    check_no_whitespace(misc_plates_d, "miscellaneous_plate names")
    for plate, plate_dict in misc_plates_d.items():
        misc_plates[plate] = {}
        if not req_keys.issubset(plate_dict):
            raise ValueError(f"miscellaneous_plate {plate} lacks {req_keys=}")
        misc_plates[plate]["viral_library"] = plate_dict["viral_library"]
        misc_plates[plate]["neut_standard_set"] = plate_dict["neut_standard_set"]
        samples = pd.read_csv(plate_dict["samples_csv"])
        if not {"well", "fastq"}.issubset(samples.columns):
            raise ValueError(
                f"{plate_dict['samples_csv']} lacks columns 'well', 'fastq'"
            )
        if len(samples) != samples["well"].nunique():
            raise ValueError(
                f"{plate_dict['samples_csv']} has non-unique entries in 'well' column"
            )
        check_no_whitespace(samples["well"], f"miscellaneous_plate {plate} wells")
        misc_plates[plate]["wells"] = samples.set_index("well")["fastq"].to_dict()

        # get default barcode parser params, update if specified per plate
        misc_plates[plate]["illumina_barcode_parser_params"] = copy.deepcopy(
            config["illumina_barcode_parser_params"]
        )
        if "illumina_barcode_parser_params" in plate_dict:
            misc_plates[plate]["illumina_barcode_parser_params"].update(
                plate_dict["illumina_barcode_parser_params"]
            )

    return misc_plates


def process_plate(plate, plate_params):
    """Process a plot from the configuration."""

    # Process plate parameters
    req_plate_params = {
        "group",
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
    check_no_whitespace(
        (
            value
            for drops in plate_params["manual_drops"].values()
            for drop in (drops or [])
            for value in (drop if isinstance(drop, list) else [drop])
        ),
        f"{plate=} 'manual_drops' values",
    )

    # a manual drop of a barcode in just some wells cannot be combined with collapsing
    # the barcodes for a strain, as the strain's collapsed counts would then be summed
    # over a different set of barcodes in those wells than in the rest of the plate
    if get_collapse_strain_barcodes(config):
        invalid_drops = {
            drop_type: drops
            for (drop_type, drops) in plate_params["manual_drops"].items()
            if drops and (drop_type in {"barcode_wells", "barcode_serum_replicates"})
        }
        if invalid_drops:
            raise ValueError(
                f"{plate=} cannot use 'collapse_strain_barcodes' with the "
                f"'manual_drops' of {sorted(invalid_drops)}, as the collapsed counts "
                "for a strain would then be summed over a different set of barcodes in "
                "those wells than in the rest of the plate, biasing its fraction "
                "infectivity. Instead use 'barcodes' to drop a barcode from the entire "
                f"plate.\n{invalid_drops=}"
            )

    plate_d = copy.deepcopy(plate_params)
    plate_d["group"] = str(plate_d["group"])
    plate_d["date"] = str(plate_d["date"])
    if not re.fullmatch(r"\d{4}\-\d{2}\-\d{2}", str(plate_d["date"])):
        raise ValueError(f"{plate=} {plate_d['date']=} not in YYYY-MM-DD format")

    # get default barcode parser params, update if specified per plate
    plate_d["illumina_barcode_parser_params"] = copy.deepcopy(
        config["illumina_barcode_parser_params"]
    )
    if "illumina_barcode_parser_params" in plate_params:
        plate_d["illumina_barcode_parser_params"].update(
            plate_params["illumina_barcode_parser_params"]
        )

    # determine whether this plate uses 'dilution_factor' or 'concentration'
    dilution_factor_or_concentration = plate_params.get(
        "dilution_factor_or_concentration", "dilution_factor"
    )
    if dilution_factor_or_concentration not in {"dilution_factor", "concentration"}:
        raise ValueError(
            f"{plate=} has invalid {dilution_factor_or_concentration=}; must be "
            "'dilution_factor' or 'concentration'"
        )
    concentration_units = plate_params.get("concentration_units")
    if dilution_factor_or_concentration == "concentration":
        if not concentration_units:
            raise ValueError(
                f"{plate=} uses 'concentration' but lacks a 'concentration_units' value"
            )
    elif concentration_units is not None:
        raise ValueError(
            f"{plate=} specifies 'concentration_units' but is not using 'concentration' "
            "(set 'dilution_factor_or_concentration: concentration' to use it)"
        )
    plate_d["dilution_factor_or_concentration"] = dilution_factor_or_concentration
    plate_d["concentration_units"] = concentration_units

    # Process samples_csv to create the sample data frame; the value of
    # 'dilution_factor_or_concentration' is also the name of the samples_csv column
    req_sample_cols = [
        "well",
        "serum",
        dilution_factor_or_concentration,
        "replicate",
        "fastq",
    ]
    samples_df = pd.read_csv(
        plate_params["samples_csv"], comment="#", dtype={"serum": str}
    )
    if {"dilution_factor", "concentration"}.issubset(samples_df.columns):
        raise ValueError(
            f"{plate=} 'samples_csv' has both 'dilution_factor' and 'concentration' "
            "columns; include only the one matching 'dilution_factor_or_concentration'"
        )
    if not set(req_sample_cols).issubset(samples_df.columns):
        raise ValueError(f"{plate=} {samples_df.columns=} lacks {req_sample_cols=}")

    if samples_df["serum"].isnull().any():
        raise ValueError(f"{plate=} 'samples_csv' has null values in 'serum' column")

    # try to turn columns of ints and NAs into Int64 to avoid ints appearing as floats
    # (only cosmetic; concentrations stay float and are not coerced)
    coerce_int_cols = ["replicate"]
    if dilution_factor_or_concentration == "dilution_factor":
        coerce_int_cols.append("dilution_factor")
    for col in coerce_int_cols:
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
                    row["serum"]
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
                        if pd.isnull(row[dilution_factor_or_concentration])
                        else f"_{row[dilution_factor_or_concentration]}"
                    )
                ),
                axis=1,
            ),
            sample=lambda x: x["sample_noplate"].map(lambda s: f"{plate}_{s}"),
            plate=plate,
            plate_replicate=lambda x: x.apply(
                lambda row: (
                    plate
                    + ("" if row["one_serum_replicate"] else f"{- row['replicate']}")
                ),
                axis=1,
            ),
        )
        .drop(columns="one_serum_replicate")
    )

    check_no_whitespace(
        samples_df["serum_replicate"], f"{plate=} serum/replicate names"
    )
    check_no_whitespace(samples_df["well"], f"{plate=} wells")

    duplicated_samples = samples_df[samples_df.duplicated("sample", False)]
    if len(duplicated_samples):
        raise ValueError(f"Duplicated samples for {plate=}:\n{duplicated_samples}")

    # make sure serum_replicate and dilution_factor / concentration are unique
    dup_rows = (
        samples_df.assign(
            duplicates=lambda x: (
                x.groupby(
                    ["serum_replicate", dilution_factor_or_concentration],
                    dropna=False,
                )["sample"].transform("count")
            ),
        )
        .query("duplicates > 1")
        .drop(columns="duplicates")
    )
    if len(dup_rows):
        raise ValueError(f"{plate=} has duplicated serum / replicates:\n{dup_rows}")

    # make sure dilution_factor / concentration is valid; no-serum wells must be
    # specified via serum == "none" (blank value), NOT by setting the value to 0
    if dilution_factor_or_concentration == "dilution_factor":
        if not (
            (samples_df["dilution_factor"] >= 1) | (samples_df["serum"] == "none")
        ).all():
            raise ValueError(
                f"{plate=} has dilution factors not >= 1 for non-'none' serum. "
                "Specify no-serum wells by setting 'serum' to 'none' (leaving "
                "'dilution_factor' blank), not by setting 'dilution_factor' to 0."
            )
    else:
        if not (
            (samples_df["concentration"] > 0) | (samples_df["serum"] == "none")
        ).all():
            raise ValueError(
                f"{plate=} has concentrations not > 0 for non-'none' serum. "
                "Specify no-serum wells by setting 'serum' to 'none' (leaving "
                "'concentration' blank), not by setting 'concentration' to 0."
            )

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


def get_plate_comparators(plates):
    """Get dict keyed by plate with values the plates its titers can be correlated with.

    A plate's titers can be correlated with those of another plate in the same group
    that shares at least one serum with it. Titers are only correlated within a group,
    as different groups can report titers in different units.

    Computed from the sample tables in the configuration, and so before any QC. QC can
    only remove sera, so a plate with no comparator here cannot gain one, while a plate
    with comparators here can lose all of its shared titers to QC (which the plate's
    report notes).

    """
    plate_sera = {
        plate: set(plate_d["samples"]["serum"]) - {"none"}
        for (plate, plate_d) in plates.items()
    }
    return {
        plate: [
            other
            for other in plates  # config order, which orders a plate's comparators
            if (other != plate)
            and (plates[other]["group"] == plate_d["group"])
            and (plate_sera[plate] & plate_sera[other])
        ]
        for (plate, plate_d) in plates.items()
    }


@functools.lru_cache
def groups_sera_plates():
    """Get dict keyed by /groupserum with values lists of plates with titers."""
    csv_file = checkpoints.groups_sera_by_plate.get().output.csv
    return (
        pd.read_csv(csv_file, dtype={"group": str, "serum": str})
        .assign(plates=lambda x: x["plates"].str.split(";"))
        .set_index(["group", "serum"])["plates"]
        .to_dict()
    )
