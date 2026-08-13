"""Make the modules in ``scripts/`` importable by the tests.

When ``snakemake`` runs a ``script:``, it adds that script's directory to ``sys.path`` so
the script can import its siblings. The tests import those modules directly, so they have
to do the same.

"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "scripts"))
