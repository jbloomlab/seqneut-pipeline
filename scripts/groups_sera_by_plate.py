"""Get list of groups/sera and plates for each."""

import sys

import pandas as pd


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

plates = snakemake.params.plates
curvefit_csvs = snakemake.input.csvs

assert len(plates) == len(curvefit_csvs) == len(set(plates))

dfs = []
for plate, curvefit_csv in zip(plates, curvefit_csvs):
    dfs.append(
        pd.read_csv(curvefit_csv)[["group", "serum"]]
        .drop_duplicates()
        .assign(plate=plate)
    )

(
    pd.concat(dfs)
    .groupby(["group", "serum"], as_index=False)
    .aggregate(plates=pd.NamedAgg("plate", lambda s: ";".join(sorted(s))))
    .to_csv(snakemake.output.csv, index=False)
)
