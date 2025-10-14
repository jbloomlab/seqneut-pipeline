# CLAUDE.md - Context for seqneut-pipeline

## Repository Overview

This is `seqneut-pipeline`, a modular Snakemake pipeline for analyzing high-throughput sequencing-based neutralization assays developed in the [Bloom lab](https://jbloomlab.org). The pipeline processes FASTQ files from barcoded viral libraries to compute neutralization titers for sera against viruses.

**Key Publications**:
- [Loes et al (2024)](https://doi.org/10.1128/jvi.00689-24) describes these assays.
- [Kikawa et al (2025a)](https://elifesciences.org/reviewed-preprints/106811) uses the assays.
- [Kikawa et al (2025b)](https://www.biorxiv.org/content/10.1101/2025.09.06.674661v1) uses the assays at scale.

**Usage Model**: This pipeline is designed to be included as a **git submodule** in project-specific repositories, not used standalone.

## Core Workflow

The pipeline performs these sequential steps:

1. **Barcode Counting** (`count_barcodes` rule)
   - Parses FASTQ files to count viral and neutralization-standard barcodes
   - Uses `dms_variants.illuminabarcodeparser.IlluminaBarcodeParser`
   - Outputs: counts, invalid barcodes, and read fates per sample

2. **Plate Processing** (`process_plate` rule)
   - Applies extensive QC thresholds to drop low-quality data
   - Computes fraction infectivity for each viral barcode at each serum concentration
   - Fits Hill-curve neutralization curves using `neutcurve` package
   - Outputs: fraction infectivity CSV, curve fits (CSV + pickle), QC drops YAML

3. **Serum-Level Aggregation** (`group_serum_titers` rule)
   - Aggregates curve fits across replicates (potentially from multiple plates)
   - Computes median titers across replicates
   - Applies QC to identify outlier replicates
   - Outputs: per-replicate titers, median titers, curves PDF, QC drops YAML

4. **Final Aggregation** (`aggregate_titers` rule)
   - Combines all sera titers for each group
   - Creates interactive visualizations
   - Outputs: aggregated titers CSV and pickled `CurveFits` objects

5. **Documentation** (`build_docs` rule)
   - Converts Jupyter notebooks to HTML
   - Builds navigable documentation for GitHub Pages
   - Outputs: `docs/` directory with index.html and all results

## Key Files and Their Roles

### Pipeline Core Files
- **`seqneut-pipeline.smk`**: Main Snakemake rules file (included by project Snakefiles)
- **`funcs.smk`**: Python functions for processing configuration (plate setup, sample validation)
- **`notebook_funcs.py`**: Helper functions used inside Jupyter notebooks (e.g., `display_curve_fig`); the main purpose of this function currently is to display very large matplotlib figures in a way where the rendered HTML of the Jupyter notebooks does not take up too much space.
- **`environment.yml`**: Conda environment specification for the pipeline

### Scripts (in `scripts/`)
- **`count_barcodes.py`**: Counts barcodes from FASTQ files
- **`groups_sera_by_plate.py`**: Creates mapping of sera to plates (checkpoint)
- **`build_docs.py`**: Builds HTML documentation with markdown-based index

### Notebooks (in `notebooks/`)
- **`process_plate.py.ipynb`**: Processes a single plate (QC, fraction infectivity, curve fitting)
- **`group_serum_titers.py.ipynb`**: Aggregates titers for a serum across plates
- **`aggregate_titers.py.ipynb`**: Combines all titers and creates visualizations
- **`aggregate_qc_drops.py.ipynb`**: Summarizes all QC filtering

### Test Example
- **`test_example/`**: Small working example with subsetted data
- **`test_example/config.yml`**: Example configuration (QC thresholds set more leniently for small data)
- **`test_example/Snakefile`**: Shows how to include and run the pipeline

## Configuration Structure

The pipeline is configured via a `config.yml` file in the parent repository. Key sections:

### Top-Level Keys

1. **`seqneut-pipeline`**: Path to this submodule (typically `seqneut-pipeline`)
2. **`docs`**: Path for HTML output (typically `docs`)
3. **`description`**: Markdown description for documentation
4. **`viral_libraries`**: Dict mapping library names to CSV files with barcode→strain mappings
5. **`neut_standard_sets`**: Dict mapping set names to CSV files with neutralization standard barcodes
6. **`illumina_barcode_parser_params`**: Global parameters for barcode parsing
7. **`curve_display_method`**: How to display curves in notebooks (`inline`, `png8`, or `no_display`)

### Plate Configuration

Each plate under `plates:` has:
- **`group`**: Categorical label (e.g., "serum" or "pilot")
- **`date`**: YYYY-MM-DD format
- **`viral_library`**: Name from `viral_libraries`
- **`neut_standard_set`**: Name from `neut_standard_sets`
- **`samples_csv`**: CSV with columns: well, serum, dilution_factor, replicate, fastq
- **`manual_drops`**: Explicit drops (wells, barcodes, barcode_wells, barcode_serum_replicates, serum_replicates)
- **`qc_thresholds`**: QC filters for plate processing (see QC section below)
- **`curvefit_params`**: Parameters for curve fitting (fixtop, fixbottom, fixslope, frac_infectivity_ceiling)
- **`curvefit_qc`**: QC for curve fits (max_frac_infectivity_at_least, goodness_of_fit)
- **`illumina_barcode_parser_params`** (optional): Plate-specific overrides

### Serum Titer Configuration

- **`default_serum_titer_as`**: How to compute titers (`midpoint` or `nt50`)
- **`default_serum_qc_thresholds`**: QC for serum-virus pairs (min_replicates, max_fold_change_from_median)
- **`sera_override_defaults`**: Per-serum/group overrides for titer calculation or QC

### Miscellaneous Plates (Optional)

- **`miscellaneous_plates`**: Plates to just count barcodes without curve fitting

## Quality Control (QC) System

The pipeline applies multi-stage QC to ensure data quality:

### Plate-Level QC (`qc_thresholds`)

1. **`avg_barcode_counts_per_well`**: Minimum average barcode counts per well
2. **`min_neut_standard_frac_per_well`**: Minimum fraction of counts from neut standard
3. **`no_serum_per_viral_barcode_filters`**: Viral barcode consistency in no-serum samples
   - `min_frac`: Minimum fraction across wells
   - `max_fold_change`: Maximum deviation from median
   - `max_wells`: Allowed failures before dropping barcode
4. **`per_neut_standard_barcode_filters`**: Neut standard barcode consistency
5. **`min_neut_standard_count_per_well`**: Minimum total neut standard counts
6. **`min_no_serum_count_per_viral_barcode_well`**: Minimum viral barcode counts in no-serum
7. **`max_frac_infectivity_per_viral_barcode_well`**: Maximum allowed fraction infectivity
8. **`min_dilutions_per_barcode_serum_replicate`**: Minimum dilution points for curve fitting

### Curve Fit QC (`curvefit_qc`)

1. **`max_frac_infectivity_at_least`**: Minimum max fraction infectivity (set to 0 to allow fully neutralized)
2. **`goodness_of_fit`**: Curve must pass EITHER R² OR RMSD threshold
   - `min_R2`: Minimum coefficient of determination
   - `max_RMSD`: Maximum root-mean-square deviation
3. **`serum_replicates_ignore_curvefit_qc`**: List of serum-replicates to exempt
4. **`barcode_serum_replicates_ignore_curvefit_qc`**: List of [barcode, serum_replicate] pairs to exempt

### Serum-Level QC (`default_serum_qc_thresholds`)

1. **`min_replicates`**: Minimum replicates for a titer
2. **`max_fold_change_from_median`**: Maximum fold-change between replicate and median
3. **`viruses_ignore_qc`**: Viruses exempt from QC for ALL sera

## Important Data Structures

### Sample Organization

Samples are uniquely identified by:
- **`sample`**: `{plate}_{serum_replicate}_{dilution_factor}` (or without dilution for no-serum)
- **`serum_replicate`**: `{serum}` or `{serum}-{replicate}` if multiple replicates exist
- **`plate_barcode`**: `{plate_replicate}-{barcode}` for tracking replicates across plates

### Fraction Infectivity Calculation

For viral barcode `v` in sample `s`:

```
F(v,s) = [c(v,s) / Σ(neut_standard_counts(s))] / median_over_no_serum_samples[c(v,s0) / Σ(neut_standard_counts(s0))]
```

This normalizes each sample by neutralization standard counts, then by the median ratio in no-serum samples.

### Curve Fitting

Uses `neutcurve.CurveFits` to fit Hill curves with configurable:
- **`fixtop`**: Top plateau (can be fixed value, range [min, max], or free)
- **`fixbottom`**: Bottom plateau (typically 0)
- **`fixslope`**: Slope (typically constrained to "reasonable" range like [0.8, 10] for serum dilutions)
- **`frac_infectivity_ceiling`**: Ceiling applied before fitting (typically 1)

Titers computed as:
- **`nt50`**: Reciprocal of IC50 (concentration at 50% infectivity)
- **`midpoint`**: Reciprocal of curve midpoint (differs from NT50 when plateaus ≠ 0 and 1)

## Output Files

### Key Results (should be tracked in git)
- **`results/aggregated_titers/titers_{group}.csv`**: Final titers (median across replicates) after all QC
- **`results/aggregated_titers/curvefits_{group}.pickle`**: Pickled `neutcurve.CurveFits` objects
- **`results/sera/{group}_{serum}/titers.csv`**: Per-serum median titers
- **`results/sera/{group}_{serum}/titers_per_replicate.csv`**: Individual replicate titers
- **`results/plates/{plate}/frac_infectivity.csv`**: Fraction infectivity per plate
- **`results/plates/{plate}/curvefits.csv`**: Curve fit parameters per plate
- **`results/qc_drops/*.yml`**: Summaries of QC filtering
- **`docs/`**: HTML documentation for GitHub Pages

### Intermediate Files (may not need git tracking)
- **`results/barcode_counts/`**: Per-sample barcode counts
- **`results/barcode_fates/`**: Read fate statistics
- **`results/plates/{plate}/*.ipynb`**: Jupyter notebooks (HTML versions in docs)
- **`results/sera/{group}_{serum}/*.ipynb`**: Jupyter notebooks (HTML versions in docs)

## Using the Pipeline in a New Project

1. Create a new repository for your project
2. Add seqneut-pipeline as a submodule:
   ```bash
   git submodule add https://github.com/jbloomlab/seqneut-pipeline
   ```
3. Create `config.yml` with your experiment configuration
4. Create top-level `Snakefile`:
   ```python
   import os
   configfile: "config.yml"
   include: os.path.join(config["seqneut-pipeline"], "seqneut-pipeline.smk")

   rule all:
       input:
           seqneut_pipeline_outputs
   ```
5. Create conda environment and run:
   ```bash
   conda env create -f seqneut-pipeline/environment.yml
   conda activate seqneut-pipeline
   snakemake -j <jobs> --software-deployment-method conda
   ```

## Real-World Example

See [flu-seqneut-2025](https://github.com/jbloomlab/flu-seqneut-2025) for a complete real-world example analyzing human neutralizing antibody responses to influenza H1N1 and H3N2 viruses.

## Common Tasks for Claude

### When Adding/Modifying Plates
1. Update `config.yml` with new plate configuration
2. Ensure `samples_csv` has correct format (well, serum, dilution_factor, replicate, fastq)
3. Check QC thresholds are appropriate for the data
4. After running, examine QC drops in `results/qc_drops/` and adjust if needed

### When Troubleshooting QC Failures
1. Look at `results/qc_drops/plate_qc_drops.yml` for plate-level drops
2. Look at `results/qc_drops/barcode_qc_drops.yml` for barcode-specific drops
3. Examine notebooks in `docs/` for interactive plots showing which data failed
4. Consider adding to `manual_drops` or adjusting QC thresholds if failures seem incorrect

### When Fitting Curves Poorly
1. Check `curvefit_params` - ensure `fixtop`, `fixbottom`, `fixslope` are reasonable
2. Look at curve plots in plate notebooks to diagnose issues
3. Consider adjusting `frac_infectivity_ceiling` if many points exceed 1
4. Use `barcode_serum_replicates_ignore_curvefit_qc` for specific problematic fits

### When Adding Custom Analysis
1. Add custom rules to top-level `Snakefile`
2. Use `results/aggregated_titers/curvefits_{group}.pickle` for making custom plots
3. Add custom HTML outputs to `add_htmls_to_docs` dict in top-level `Snakefile`

## Code Style and Testing

- Python code formatted with [black](https://github.com/psf/black)
- Snakemake code formatted with [snakefmt](https://github.com/snakemake/snakefmt)
- Linted with [ruff](https://github.com/astral-sh/ruff) and `snakemake --lint`
- Tested via GitHub Actions on `test_example/`

## Important Notes

1. **Groups cannot contain underscores or pipes** (causes wildcard issues)
2. **Dates must be in YYYY-MM-DD format**
3. **Barcodes must all be same length** within a library
4. **Each plate must have at least one no-serum sample** (serum = "none")
5. **QC thresholds in test_example are lenient** - real experiments should use stricter values
6. **Always examine QC drops** after adding new data before trusting results
7. **Curve fitting uses data WITH ceiling applied** but raw values are also saved
