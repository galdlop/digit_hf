# Reanalysis of the DIGIT-HF Trial: Time-Varying Treatment Effects and Violation of Proportional Hazards

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Made with R](https://img.shields.io/badge/Made%20with-R-276DC3.svg)]()
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

This repository contains the **exact reproducible analysis** used in the medRxiv manuscript:

**_Time-Varying Treatment Effects and Violation of Proportional Hazards in the DIGIT-HF Trial_**  
_Preprint on medRxiv (2025)._  
👉 *(Link to be added once available)*

The aim is to evaluate whether the treatment effect of digitoxin in DIGIT-HF is **constant over time**, and to provide **time-varying estimands** as recommended by contemporary guidelines for survival analysis and regulatory reporting.

---

## 🔍 Overview

This repository includes:

### **1. Reconstruction of Individual Patient Data (IPD)**
- Using the validated **Guyot et al. (2012)** algorithm.  
- Reproduces the Kaplan–Meier curves from the published DIGIT-HF trial.  
- Successfully replicates the primary Cox hazard ratio (HR 0.82).

### **2. Proportional Hazards Diagnostics**
- Cox PH model  
- Schoenfeld residuals + **cox.zph** global test  
- Scaled Schoenfeld residuals vs. time  
- Log–log survival curves  
- Visual and statistical confirmation of **non-proportional hazards**

### **3. Time-Varying Hazard Ratio**
- Flexible parametric survival model (**Royston–Parmar**, `stpm2`)  
- HR(t) curve showing early benefit followed by attenuation  
- Comparison with the Cox average hazard ratio (AHR)

### **4. Landmark Estimands**
- Period-specific hazard ratios at 6, 12, 18, 24 months  
- Estimates aligned with modern estimand frameworks  
- Demonstrates time-limited benefit

### **5. Figures and Scripts**
All scripts produce the figures used in the manuscript, stored in:

---

## 📁 Repository Structure
.
├── R/
│ ├── 01_cox_PH_diagnostics.R
│ ├── 02_time_varying_HR_stpm2.R
│ ├── 03_landmark_estimands.R
│ └── utils_plotting.R
│
├── data/
│ └── reconstructed_ipd.csv
│
├── figures/
│ └── (auto-generated plots)
│
├── README.md
├── LICENSE
└── zenodo.json


---

## 🔬 Data Information (Important)

The file:



/data/reconstructed_ipd.csv


contains **reconstructed** patient-level data derived from digitized Kaplan–Meier curves using:

> Guyot P, Ades AE, Ouwens MJNM, Welton NJ.  
> *Enhanced secondary analysis of survival data: reconstructing IPD from published Kaplan–Meier curves.*  
> BMC Med Res Methodol. 2012.

These are **not original trial data** from DIGIT-HF.  
Reconstructed datasets typically achieve excellent accuracy for reproducing trial results, but small deviations in event time alignment may occur.

All analyses and figures should be interpreted accordingly.

---

## ▶️ How to Reproduce the Analysis

1. Place your reconstructed IPD file in:



/data/reconstructed_ipd.csv


with the following columns:

| Column | Description |
|--------|-------------|
| `id` | Unique subject identifier |
| `time` | Time to event or censoring (same units as trial KM curves) |
| `status` | 1 = event, 0 = censored |
| `arm` | 0 = control/placebo, 1 = digitoxin |

2. Run any of the scripts in `/R/`.

The figures and diagnostic outputs will be generated automatically in `/figures/`.

---

## 📘 Citation

If you use this code or analysis in your work, please cite:

**Preprint (medRxiv):**  
_“Time-Varying Treatment Effects and Violation of Proportional Hazards in the DIGIT-HF Trial”_  
_Link to be added_

**Software / Code (Zenodo):**  


[Your Name]. Reanalysis of the DIGIT-HF Trial: Time-Varying Treatment Effects. Zenodo. DOI:10.5281/zenodo.XXXXXXX


---

## 📄 License  
This project is released under the **MIT License**, permitting reuse and adaptation with attribution.

---

## 🙌 Acknowledgements  
Thanks to the developers of `survival`, `survminer`, `rstpm2`, and all contributors to open-source R tools for survival analysis.

---

## 💬 Contact  
For questions or collaboration opportunities:  
**Your name — Your institution**  
📧 **your.email@institution.org**  
ORCID: https://orcid.org/0000-0000-0000-0000


