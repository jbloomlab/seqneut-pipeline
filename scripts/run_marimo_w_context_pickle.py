"""Run a marimo notebook with snakemake context passed via pickle."""

import os
import pickle
import subprocess
import sys

# Redirect stderr and stdout to log
log_file = open(snakemake.log[0], "w")
sys.stderr = sys.stdout = log_file

# Build context dictionary from snakemake object
context = {
    "workdir": os.getcwd(),
    "input": dict(snakemake.input),
    "output": dict(snakemake.output),
    "params": snakemake.params,
    "threads": snakemake.threads,
    "resources": dict(snakemake.resources),
    "wildcards": dict(snakemake.wildcards),
}

marimo_nb = snakemake.input.marimo_nb
marimo_html = snakemake.output.marimo_html
context_pickle = snakemake.output.context_pickle

print(f"Running marimo notebook: {marimo_nb=}")
print(f"Using context pickle: {context_pickle=}")
print(f"Context keys: {list(context.keys())=}")
for key in context:
    if context[key]:
        print(f"  For {key}, context is: {context[key]}")
print(f"Output HTML: {marimo_html=}")

# Write context to pickle file
with open(context_pickle, "wb") as f:
    pickle.dump(context, f, protocol=pickle.HIGHEST_PROTOCOL)

# Run marimo export
cmd = [
    "marimo",
    "export",
    "html",
    marimo_nb,
    "-o",
    marimo_html,
    "--",
    "--context-pickle",
    context_pickle,
]

print(f"Running command: {' '.join(cmd)}")
print("=" * 80)

try:
    # Stream output directly to log file (real-time logging)
    result = subprocess.run(
        cmd,
        stdout=log_file,
        stderr=subprocess.STDOUT,
        check=False,
    )
    print("=" * 80)
    if result.returncode != 0:
        print(f"ERROR: marimo export failed with return code {result.returncode}")
        raise subprocess.CalledProcessError(result.returncode, cmd)
    print("marimo export completed successfully")
except Exception as e:
    print(f"ERROR: Exception running marimo export: {e}")
    raise

print(f"HTML output written to: {marimo_html}")

log_file.close()
