# Test example for `seqneut-pipeline`

This subdirectory has a small example that illustrates use and testing of the `seqneut-pipeline`.

Note that this is a small example designed to run quickly, so the datasets are smaller than real ones and some thresholds have been relaxed to more lenient values in the YAML configuration than should be done for real experiments.

Also, note that the FASTQ files for the test example are contained in [./fastqs/](fastqs) to make it self contained---for real experiments, you would **not** track the FASTQ files in the repo.

The configuration for the test example is in [config.yml](config.yml), but if using this for real experiments note the comments for QC thresholds that are set to more lenient values than is ideal for real experiments.

The input data for the test example are in [./data/](data).

To run the example, build the `seqneut-pipeline` `conda` environment in [../environment.yml](../environment.yml), activate it, and then run:

    snakemake -j <n_jobs> --use-conda --keep-going

This will create the results in [./results/](results).
