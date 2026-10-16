# REMIT
Re-thinking Health Policy in Multiple Long Term Conditions

Code for the patient-level simulation in:

> Owen RK, Hazell I, Lyons J, Latimer N, Adler AI, Guthrie B, Lyons RA, Abrams KR. *SGLT-2 inhibitors across competing pathways of type 2 diabetes, chronic kidney disease and heart failure.* [Research Square, 2026, DOI – to be added on publication]

The simulation combines transition hazards from a multi-state natural history model, fitted to linked electronic health records for 1.3 million people in Wales (Wales Multimorbidity e-Cohort, SAIL Databank), with treatment effects from randomised controlled trials. It compares SGLT-2 inhibitors with standard of care across type 2 diabetes, chronic kidney disease and heart failure.

## Files

| File | Description |
| --- | --- |
| [Patient_level_simulation.R](Patient_level_simulation.R) | Patient-level simulation of 1,000,000 individuals per arm over 100 Monte Carlo iterations. Produces the transition probabilities, ages at transition and life expectancy reported in the paper. |
| [MC_stability.R](MC_stability.R) | Monte Carlo stability assessment (running-quantile trace plots and Monte Carlo standard errors), using the output of the simulation, and reported in the supplementary files. |

## Data

No individual-level data are included or needed. The model parameters, estimated within the SAIL Databank Trusted Research Environment, are written into the simulation script. Access to the underlying SAIL data is available to researchers on approval: see [Apply to work with the data](https://saildatabank.com/data/apply-to-work-with-the-data/).

## Requirements

R (version 4.5.0 used) with the packages `dplyr`, `tidyr`, `ggplot2`, `scales`, `survival` and `flexsurv`:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "scales", "survival", "flexsurv"))
```

## How to run

1. Run `Patient_level_simulation.R`. It writes:
   - `pooled_overall_incl_eligible.csv` – overall life expectancy differences
   - `pooled_transition_table.csv` – results by transition and arm
   - `pooled_difference_table.csv` – differences between arms by transition
   - `per_iteration_differences.csv` – paired differences for each iteration
2. Run `MC_stability.R` in the same working directory. It reads the first and last of these files and writes `MC_stability_summary.csv` and the stability plots.

Random seeds are set in both scripts, so results are reproducible.

## Citation

If you use this code, please cite this paper: Owen RK, Hazell I, Lyons J, Latimer N, Adler AI, Guthrie B, Lyons RA, Abrams KR. *SGLT-2 inhibitors across competing pathways of type 2 diabetes, chronic kidney disease and heart failure.* [Research Square, 2026, DOI – to be added on publication]

## License

[LICENSE](LICENSE)

## Funding

This work was supported by a Health and Care Research Wales – NIHR Advanced Fellowship (HCRW NIHR FS(A)-2023b-RO/NIHR303628) and Health Data Research UK (HDRUK2023.0030).

## Contact

Professor Rhiannon K Owen, Population Data Science, Swansea University Medical School – r.k.owen@swansea.ac.uk
