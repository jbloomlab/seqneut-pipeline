"""Tests that ``test_example/config.yml`` restates the optional keys' actual defaults.

The example configuration spells out the thresholds of
``analyze_single_well_infections`` and ``analyze_pools`` at the values they default to,
so that someone copying it can see which keys are theirs to tune. That makes the example
a second copy of values ``funcs.smk`` already holds, and a copy that nothing reads: an
entry giving a key its default is indistinguishable from one that leaves the key out. So
a default changed in ``funcs.smk`` would leave the example quietly showing the old value,
looking like a project that had set it deliberately. These tests are what makes the two
agree.

The defaults are read back out of ``process_single_well_infections`` and
``process_analyze_pools`` rather than off the ``defaults`` dicts they apply, so that what
is compared is the value an entry actually ends up with.

"""

import pathlib

import pytest
import yaml

REPO = pathlib.Path(__file__).parent.parent

# The keys the example is expected to restate: the thresholds that decide which wells are
# trusted, which dilutions are read, and what is flagged. The other optional keys are per-
# experiment facts rather than thresholds -- which wells hold the pool, what it was made
# from, which columns to show -- and the example sets those only where a pool needs them.
SINGLE_WELL_THRESHOLDS = [
    "max_other_library_frac",
    "max_unmatched_frac",
    "min_parsed_reads_per_well",
    "min_neut_standard_reads_for_titer",
]
POOL_THRESHOLDS = [
    "min_avg_barcode_count_per_well",
    "min_neut_standard_frac_per_well",
    "linearity_window_dilutions",
    "linearity_slope_tol",
    "max_stock_share_ratio",
]


@pytest.fixture(scope="module")
def example_config():
    """The example configuration, with its ``addtl_config_files`` merged in.

    The plates are in those files rather than the main one, and the sections under test
    refer to them, so the configuration is only complete once they are loaded. They hold
    separate top-level keys, which is what the pipeline requires of them, so merging is a
    plain update.

    """
    example = REPO / "test_example"
    with open(example / "config.yml") as f:
        config = yaml.safe_load(f)
    for addtl in config["addtl_config_files"]:
        with open(example / addtl) as f:
            config.update(yaml.safe_load(f))
    return config


@pytest.fixture(scope="module")
def funcs_smk(example_config):
    """``funcs.smk``'s namespace, with the globals its functions read.

    It is a Snakemake include rather than a module, so it is compiled and run here as
    the ``Snakefile`` including it would run it. Its own imports are standard library;
    ``pd`` and ``config`` are the globals it documents taking from the files that include
    it.

    """
    import pandas

    namespace = {"pd": pandas, "config": example_config}
    source = (REPO / "funcs.smk").read_text()
    # the repository's own file, run the way a `Snakefile` including it runs it:
    # there is no import that reaches a Snakemake include
    exec(compile(source, str(REPO / "funcs.smk"), "exec"), namespace)  # noqa: S102
    return namespace


@pytest.fixture(scope="module")
def misc_plates(example_config):
    """The example's plates, in the shape the ``process_`` functions take them in.

    Built here rather than through ``process_miscellaneous_plates``, which reads each
    plate's samples CSV for its wells; only the library, the standard set and the date are
    read from a plate by the functions under test.

    """
    return {
        plate: {
            "viral_library": d["viral_library"],
            "neut_standard_set": d["neut_standard_set"],
            "date": str(d["date"]) if d.get("date") else None,
        }
        for (plate, d) in example_config["miscellaneous_plates"].items()
    }


def defaults_of(process, section, misc_plates, required):
    """What `process` fills an entry's optional keys in with.

    The entry gives only `required`, so every other key of the result is that key's
    default. Taken from the function that applies the defaults rather than from the dict
    it applies, so that a default moved or overridden is still read correctly.

    """
    name, entry = next(iter(section.items()))
    minimal = {key: entry[key] for key in required}
    return process({name: minimal}, misc_plates)[name]


@pytest.mark.parametrize(
    "section_key, process_name, required, thresholds",
    [
        (
            "analyze_single_well_infections",
            "process_single_well_infections",
            ["miscellaneous_plates", "description", "analyze_by"],
            SINGLE_WELL_THRESHOLDS,
        ),
        (
            "analyze_pools",
            "process_analyze_pools",
            ["analyze_by", "miscellaneous_plates", "description"],
            POOL_THRESHOLDS,
        ),
    ],
)
def test_example_restates_defaults(
    example_config,
    funcs_smk,
    misc_plates,
    section_key,
    process_name,
    required,
    thresholds,
):
    """Every threshold the example states equals the default it is meant to show."""
    section = example_config[section_key]
    defaults = defaults_of(funcs_smk[process_name], section, misc_plates, required)
    assert set(thresholds).issubset(defaults), (
        f"{sorted(set(thresholds) - set(defaults))} are not optional keys of "
        f"{section_key}, so the example cannot be showing their defaults"
    )
    for name, entry in section.items():
        for key in thresholds:
            assert key in entry, (
                f"{section_key} {name!r} does not state {key!r}; the example states the "
                "thresholds at their defaults so that they are visible to copy"
            )
            assert entry[key] == defaults[key], (
                f"{section_key} {name!r} states {key}: {entry[key]!r}, but the default "
                f"is {defaults[key]!r}; the example shows the defaults, so change it to "
                "match or set the key deliberately somewhere the test does not check"
            )


def test_pool_linearity_defaults_match_pool_funcs(
    funcs_smk, example_config, misc_plates
):
    """``funcs.smk``'s linearity defaults are the ones ``pool_funcs`` fits with.

    The configuration's defaults are what a report is run with, and `fit_dilution_windows`
    applies its own to a caller that passes none -- the tests of it among them. Two
    different values would make the tested behavior not the configured behavior.

    """
    import pool_funcs

    defaults = defaults_of(
        funcs_smk["process_analyze_pools"],
        example_config["analyze_pools"],
        misc_plates,
        ["analyze_by", "miscellaneous_plates", "description"],
    )
    for config_key, constant in [
        ("linearity_window_dilutions", "DEFAULT_LINEARITY_WINDOW_DILUTIONS"),
        ("linearity_slope_tol", "DEFAULT_LINEARITY_SLOPE_TOL"),
    ]:
        assert defaults[config_key] == getattr(pool_funcs, constant), (
            f"the {config_key!r} default is {defaults[config_key]!r} but "
            f"`pool_funcs.{constant}` is {getattr(pool_funcs, constant)!r}"
        )
