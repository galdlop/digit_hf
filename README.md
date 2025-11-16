#Reproducible Analysis for: “Time-varying Treatment Effects and Violation of Proportional Hazards in the DIGIT-HF Trial”

This repository contains the full code, figures, and reconstructed individual-level data used to perform the statistical reanalysis reported in the following manuscript:

Time-varying Treatment Effects and Violation of Proportional Hazards in the DIGIT-HF Trial
Preprint posted on medRxiv (2025).
[Link will be added here once available]

#Overview

The objective of this repository is to provide fully transparent and reproducible code for evaluating whether the treatment effect reported in the DIGIT-HF trial is constant over time, and for illustrating how time-varying estimands can complement or refine the original Cox hazard ratio.

The analyses include:

1. Reproduction of trial curves

Reconstruction of individual patient data (IPD)
from published Kaplan–Meier curves using the Guyot et al. algorithm.

Validation by replicating the original Cox result.

2. Proportional hazards assessment

Cox proportional hazards model

Schoenfeld residuals and cox.zph PH test

Scaled Schoenfeld residuals vs. time

Log–log survival curves

3. Time-varying treatment effect

Flexible parametric survival model (Royston–Parmar, stpm2)

Estimated time-varying hazard ratio HR(t)

Comparison with overall Cox average hazard ratio (AHR)

4. Landmark estimands

Period-specific hazard ratios at 6, 12, 18, and 24 months

Interpretation consistent with modern regulatory and estimand frameworks

5. Code for all figures

All plots used in the manuscript (PH diagnostics, HR(t), log–log curves)

Exported to /figures/ during execution

Important Note on the Data

The file:

/data/reconstructed_ipd.csv


contains reconstructed individual-level survival data generated from digitized Kaplan–Meier curves, following the validated method of:

Guyot P, Ades AE, Ouwens MJNM, Welton NJ.
Enhanced secondary analysis of survival data: reconstructing IPD from published Kaplan–Meier curves.
BMC Med Res Methodol. 2012.

These are not the original patient-level data from DIGIT-HF.
Reconstructed IPD allow highly accurate replication of published results but may differ slightly in event times or censoring patterns.

All analyses should be interpreted with this in mind.
