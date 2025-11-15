# Reproducible analysis: Time-dependent treatment effect in DIGIT-HF trial 
This repository contains code and reconstructed individual-level data to reproduce the Kaplan–Meier curves reported in the DIGIT-HF trial (Bavendiek et al. in N Engl J Med 2025;393:1155-1165) and to test the proportional-hazards assumption. Analyses include:

- Reproduction of published KM curves (from digitized points).
- Cox proportional-hazards model and Schoenfeld residuals test (cox.zph).
- Scaled Schoenfeld Residuals vs. Time.
- Log-Log Survival Curves.
- Time-varying hazard ratio estimation (Cox with time interaction).
- Modeling and plotting a Landmark Analysis at 12 months.
- Plots and scripts to reproduce figures.

**Important note:** Data in `/data/reconstructed_ipd.csv` are reconstructed from published curves (Guyot et al. approach). They are not original patient-level data and should be interpreted accordingly.

## How to reproduce

1. Place your reconstructed IPD (CSV) in `/data/reconstructed_ipd.csv` with columns:
   - `id` (unique identifier)
   - `time` (time to event or censoring, same units as figure)
   - `status` (1=event, 0=censored)
   - `arm` (0=control/placebo, 1=treatment/digitoxin)
2. Run the code.
