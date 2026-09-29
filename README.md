# TFMseg

This repository contains the R code accompanying our work on change point detection and mode-identification in tensor factor models.

**Paper:** *Detection and Mode-Identification of Multiple Change Points in Tensor Factor Models*  
**arXiv:** https://arxiv.org/abs/2604.11300

## Overview

The repository provides code for:

- change point detection using TFMseg;
- mode-identification for detected change points;
- mode-informed loading space estimation;
- simulation experiments.

The complete simulation results are collected in `simulation.Rmd`. A worked example illustrating the use of TFMseg is provided in `example.Rmd`.

## Repository structure

```text
TFMseg/
├── README.md
├── simulation.Rmd
├── example.Rmd
├── data_example.rds
├── Method_funcs/
│   ├── SBS.R
│   ├── global_est.R
│   ├── TFMseg.R
│   ├── TFMseg_local_lrv.R
│   ├── mode_identification.R
│   └── packages.R
└── Simulation/
    ├── DGP.R
    ├── simu_detection.R
    ├── multi_hist.R
    ├── simu_identification.R
    ├── simu_reest.R
    └── compare_detect_LR.R
```

The main files are:

- `SBS.R`: constructs the seeded intervals used by the change point detection procedure.
- `global_est.R`: implements global PCA estimation for vector, matrix and tensor factor models, including the projected estimator.
- `TFMseg_local_lrv.R`: implements the main TFMseg change point detection procedure using interval-wise HAC/long-run variance standardisation.
- `TFMseg.R`: provides an alternative implementation using global HAC/long-run variance standardisation.
- `mode_identification.R`: implements mode-identification and the corresponding evaluation functions.
- `DGP.R`: contains the data-generating mechanism for vector, matrix and tensor factor models with multiple change points.
- `simulation.Rmd`: reproduces the complete simulation analyses and figures reported in the paper.
- `example.Rmd`: provides a worked example illustrating how to apply TFMseg to a simulated tensor time series.
- `data_example.rds`: contains the simulated dataset used in `example.Rmd`, generated under setting (S1) of the paper.

`simulation.Rmd` also uses simulation wrappers and plotting functions from the `Simulation/` directory. These helper files should therefore be kept in the repository together with the files listed above.

## Requirements

The code is written in R. The project uses the following R packages:

```r
tensorMiss
rTensor
RTFA
ggplot2
dplyr
tidyr
kableExtra
lubridate
knitr
scales
```

To render the R Markdown files, the `rmarkdown` package is also required.

Install any missing CRAN packages in the usual way, for example:

```r
install.packages(c(
  "rTensor", "ggplot2", "dplyr", "tidyr",
  "kableExtra", "lubridate", "knitr",
  "scales", "rmarkdown"
))
```

## Reproducing the simulation results

Clone or download the repository and set the repository root as the working directory.

The required functions can be sourced using repository-relative paths:

```r
dirs <- c(
  "Method_funcs",
  "Simulation"
)

for (d in dirs) {
  r_files <- list.files(d, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}
```

Then render the simulation file from the repository root:

```r
rmarkdown::render("simulation.Rmd")
```

Alternatively, open `simulation.Rmd` in RStudio and select **Knit**.

The simulation file reproduces the complete numerical experiments reported in Appendix D, *Complete numerical experiments*, of the paper.

The full simulation grid covers the settings considered in the appendix, including different change point scenarios, temporal dependence structures, change point spacings, tensor dimensions and sample sizes, as well as experiments on mode-identification, mode-informed loading space estimation and missing observations. Running the complete set of experiments may require substantial computation.

## Example

For a worked example illustrating change point detection and mode-identification using TFMseg, see `example.Rmd`. The example uses the simulated dataset stored in `data_example.rds`.
