# Fit BD Model to Truncated PBD-Simulated Trees with Trimmed Terminal Branches

## Overview

This pipeline fits a Birth-Death (BD) model to truncated phylogenies based on trees simulated with the Protracted Birth-Death (PBD) model. The scripts are designed to run multiple replicates and combinations of parameters, primarily intended to be executed on a computing cluster.

## Files

* `1-fit_truncated_PBD_trees.Rmd`: R Markdown file to fit a BD model adapted to truncated phylogenies on PBD-simulated trees.
* `1-fit_truncated_PBD_trees.R`: Script version of the `.Rmd` file.
* `1-fit_truncated_PBD_trees.nb.html`: HTML report generated from the `1-fit_truncated_PBD_trees.Rmd` file.
* `1-submission.sub`: Job submission script for running the fitting process on a cluster.
* `2-combine_CSVs.R`: Script to combine the output CSV files (stored in "simulations_output/4-Truncated_PBD") from multiple replicates into a single CSV file for further analysis.
* `3-plot_comparison.R`: Script to generate comparison plots from the combined CSV data. These plots are used to visualize and compare the results of the model fitting process.

## Usage

### Fitting the BD Model

To fit the BD model to truncated PBD-simulated trees, run the following script on your local machine:
```
Rscript 1-fit_truncated_PBD_trees.R [parameter_row]
```
The *parameter_row* argument corresponds to the row in the parameter file "simulations_output/1-PBD/all_simulations_inference.csv".

### Combining CSV Files
After running multiple replicates, combine the output CSV files into a single file using:

```
Rscript 2-combine_CSVs.R
```
The result will be stored as "simulations_output/4-Truncated_PBD/combined_long_df.csv".

### Generating Comparison Plots

To generate comparison plots from the combined CSV data, run:

```
Rscript 3-plot_comparison.R
```
The resulting figures will be stored as "fig/SM_truncation_lambda.pdf" and "fig/SM_truncation_mu.pdf".
