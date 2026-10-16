# Drought-to-Flood Event Detection for Pitfalls HESS paper: https://hess.copernicus.org/articles/29/6069/2025/

This repository contains the R code used to identify droughts, floods, and drought-to-flood transitions in daily streamflow records. The scripts are provided to support reproducibility of the event-detection methodology and associated analyses.

## Repository structure

```text
├── code/
│   └── multi_functions_copy.R
├── case_study_data/
│   └── combined_flow_data.csv
├── reproduce_events.R
└── README.md
```

- `reproduce_events.R` – Main script for event detection, daily classification, and visualization.
- `code/multi_functions_copy.R` – Supporting functions for event identification.
- `case_study_data/combined_flow_data.csv` – Input daily streamflow data. This is for demonstration only. I do not own the original data. 

## Requirements

The analysis was implemented in R and requires the following packages:

```r
install.packages(c(
  "dplyr", "data.table", "evd", "zoo", "ggplot2"
))
```

## Usage

1. Download or clone the repository.
2. Set `project_dir` in `reproduce_events.R` to the repository's root directory.
3. Run `reproduce_events.R`.

The input CSV must contain the columns `ID` (catchment identifier), `Date` (YYYY-MM-DD), and `streamflow` (daily streamflow).

## Methodology

The script applies six combinations of fixed, variable, and seasonal drought and flood thresholds. The analysis uses the original method-specific percentiles, a smoothing parameter of 31 observations, a ±15-day variable-threshold window, a minimum event length of 30 days, and a minimum recovery period of 15 days.

## Outputs

Results are saved to `reproduction_results/` and include:

- Daily streamflow classifications and associated thresholds (CSV).
- Identified event indices (CSV).
- Catchment-level time-series visualizations (JPEG).

The supporting event-detection functions are used without modification relative to what was used for the original manuscript.

## Citation

If you use this code, please cite the associated publication https://hess.copernicus.org/articles/29/6069/2025/.
