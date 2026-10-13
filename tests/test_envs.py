"""Check that the ``conda`` environments agree on the version of every package they share.

Each rule runs in the minimal environment for its script, defined in ``envs/``, so that
changing one dependency rebuilds only the environments that use it. The top-level
``environment.yml`` is the environment the pipeline is run and developed from, and holds
everything the tests import. A package listed in more than one of these files must be
pinned the same way in each, or a script could pass its tests against one version and run
in the pipeline against another. This guards only what the files list: a package left
unpinned, or pulled in only as a dependency, is solved separately in each environment and
can differ between them.

"""

import collections
import pathlib
import re

from ruamel.yaml import YAML

REPO = pathlib.Path(__file__).parent.parent

ENV_DIR = REPO / "envs"

ENV_FILES = [
    REPO / "environment.yml",
    *sorted([*ENV_DIR.glob("*.yml"), *ENV_DIR.glob("*.yaml")]),
]


def specs(path):
    """Map each package in the environment file at `path` to its spec, as written."""
    deps = YAML(typ="safe").load(path)["dependencies"]
    flat = []
    for dep in deps:
        if isinstance(dep, dict):
            flat.extend(dep["pip"])
        else:
            flat.append(dep)
    # A channel prefix (`conda-forge::pandas`) is not part of the name, and `pip` treats
    # `-` and `_` in a name as the same, so `dms-variants` is `dms_variants`
    names = [
        re.split(r"[=<>!~ ]", dep.split("::")[-1], maxsplit=1)[0]
        .lower()
        .replace("-", "_")
        for dep in flat
    ]
    assert len(set(names)) == len(names), f"{path.name} lists a package twice: {flat}"
    return dict(zip(names, flat, strict=True))


def test_found_env_files():
    """The glob found the per-rule environments, so the other test is not vacuous."""
    assert len(ENV_FILES) > 1


def test_shared_packages_are_pinned_identically():
    """A package in several environment files has the same spec in each."""
    by_package = collections.defaultdict(dict)
    for path in ENV_FILES:
        for package, spec in specs(path).items():
            by_package[package][path.relative_to(REPO).as_posix()] = spec
    mismatched = {
        package: files
        for package, files in by_package.items()
        if len(set(files.values())) > 1
    }
    assert not mismatched, mismatched


def test_dev_environment_holds_every_package():
    """`environment.yml` holds every package of the per-rule environments.

    The tests import the scripts' modules directly, so they need every package any rule
    does, and `environment.yml` is the environment they run in.

    """
    dev = specs(REPO / "environment.yml")
    for path in ENV_FILES[1:]:
        missing = specs(path).keys() - dev.keys()
        assert not missing, f"{path.name}: {missing}"


def test_rules_use_exactly_the_files_in_envs():
    """Every `conda:` path in `seqneut-pipeline.smk` exists, and every env file is used.

    `snakemake -n` checks neither without `--software-deployment-method conda`, so a
    misspelled path otherwise surfaces only when the pipeline is run.

    """
    smk = (REPO / "seqneut-pipeline.smk").read_text()
    used = set(re.findall(r'conda:\s*"([^"]+)"', smk))
    assert used == {path.relative_to(REPO).as_posix() for path in ENV_FILES[1:]}
