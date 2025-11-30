# Revisiting the DIGIT-HF Trial: A Methodological Re-analysis Demonstrating Time-Dependent Treatment Effects

[![License: MIT](https://img.shields.io/badge/Code-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![License: CC BY-NC 4.0](https://licensebuttons.net/l/by/4.0/88x31.png)](https://creativecommons.org/licenses/by-nc/4.0/deed.en)

## 📋 Overview

This repository contains the **complete reproducible analysis** for our methodological re-analysis of the DIGIT-HF trial, available as a preprint on medRxiv: https://www.medrxiv.org/content/10.1101/2025.11.22.25339857v1
.

**Citation:**
```
Aldama-López G, López-Vázquez D, Rebollal-Leal F. Revisiting the DIGIT-HF Trial: 
A Methodological Re-analysis Demonstrating Time-Dependent Treatment Effects. 
medRxiv. 2025. doi: 10.1101/2025.11.22.25339857
```

**Original trial:**
```
Bavendiek U, et al. Digitoxin in Patients with Heart Failure and Reduced Ejection 
Fraction. N Engl J Med. 2025;393(12):1155-1165. doi: 10.1056/NEJMoa2408504
```

## 🎯 Key Findings

Our re-analysis reveals that the DIGIT-HF trial violates the proportional hazards assumption (p=0.019), demonstrating that:

- ✅ **Early benefit** (0-18 months): HR 0.69 (95% CI: 0.54-0.88), p<0.001
- ❌ **No late benefit** (>18 months): HR 0.99 (95% CI: 0.77-1.28), p=0.94
- ⚠️ The reported constant HR of 0.82 masks this temporal heterogeneity

## 📊 Repository Contents

### Data
- `/data/reconstructed_ipd.csv` - Reconstructed individual patient data using the Guyot et al. (2012) method
  - **Important:** These are reconstructed data from published Kaplan-Meier curves, not original trial data

### Code
- `digit_hf_git_hub.R` - Complete analysis pipeline including:
  - Reproduction of published Kaplan-Meier curves
  - Cox proportional hazards model
  - Formal testing of PH assumption (Schoenfeld residuals)
  - Log-log survival curves
  - Time-dependent hazard ratio modeling (Royston-Parmar flexible parametric models)
  - Landmark analysis at 18 months

### Figures
All figures can be reproduced by running the R script:
- **Figure 1:** Replication of cumulative incidence curves
- **Figure 2:** Log-log survival curves
- **Figure 3:** Scaled Schoenfeld residuals
- **Figure 4:** Time-dependent hazard ratio
- **Figure 5:** Landmark analysis (early vs late periods)
- **Figure 6:** Graphical abstract (combined)

## 🔧 Requirements

### R Version
- R ≥ 4.0.0

### Required Packages
```r
# Install all required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse,    # Data manipulation and plotting (≥ 2.0.0)
  survival,     # Survival analysis (≥ 3.5-0)
  survminer,    # Survival plots (≥ 0.4.9)
  rstpm2,       # Flexible parametric survival models (≥ 1.6.0)
  broom,        # Tidy model outputs (≥ 1.0.0)
  patchwork,    # Combine plots (≥ 1.1.0)
  scales,       # Plot scales (≥ 1.2.0)
  showtext      # Custom fonts (optional)
)
```

## 🚀 How to Reproduce

### Quick Start
```r
# Clone the repository
git clone https://github.com/galdlop/digit_hf.git
cd digit_hf

# Open R/RStudio and run
source("digit_hf_git_hub.R")
```

### Step-by-Step

1. **Clone or download this repository**
```bash
   git clone https://github.com/galdlop/digit_hf.git
```

2. **Ensure data file is in place**
   - The reconstructed IPD should be at: `data/reconstructed_ipd.csv`
   - Columns required:`time`, `event`, `arm`

3. **Run the analysis**
```r
   # Open R or RStudio
   # Set working directory to the repository folder
   setwd("path/to/digit_hf")
   
   # Run the complete analysis
   source("digit_hf_git_hub.R")
```

4. **Output**
   - All figures will be generated and displayed
   - High-resolution PNG files will be saved automatically
   - Statistical test results will be printed to console

## 📈 Data Structure

The reconstructed IPD file (`reconstructed_ipd.csv`) contains the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `time` | numeric | Time to event or censoring (in months) |
| `event` | integer | Event indicator (1=event occurred, 0=censored) |
| `arm` | integer | Treatment arm (0=placebo, 1=digitoxin) |


## 📝 Methods

### Data Reconstruction
Individual patient data were reconstructed from published Kaplan-Meier curves using the algorithm described by:
> Guyot P, et al. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol. 2012;12:9.

### Statistical Analysis
- **Cox proportional hazards model:** Standard time-to-event analysis
- **Schoenfeld residuals test:** Formal test for PH assumption (cox.zph function)
- **Flexible parametric models:** Royston-Parmar models with time-varying coefficients (rstpm2 package)
- **Landmark analysis:** Conditional analysis at 18 months to separate early vs late effects

All analyses were conducted in R version 4.2. Complete session info available in the script output.

## ⚠️ Important Limitations

1. **Reconstructed data:** The IPD used in this analysis are reconstructed from published curves, not original trial data
2. **Limited variables:** Only time, event status, and treatment arm are available
3. **No access to:** Baseline characteristics, biomarkers, adherence data, or competing risks information
4. **Landmark choice:** The choice of 18 months as the landmark time was pre-specified based on the time-dependent HR crossing the null effect

## 📚 References

1. Bavendiek U, et al. Digitoxin in Patients with Heart Failure and Reduced Ejection Fraction. N Engl J Med. 2025;393(12):1155-1165.

2. Guyot P, et al. Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves. BMC Med Res Methodol. 2012;12:9.

3. Cox DR. Regression Models and Life-Tables. J R Stat Soc Series B. 1972;34(2):187-220.

4. Royston P, Parmar MK. Flexible parametric proportional-hazards and proportional-odds models for censored survival data. Stat Med. 2002;21(15):2175-97.

5. Stensrud MJ, Hernán MA. Why Test for Proportional Hazards? JAMA. 2020;323(14):1401-1402.

## 👥 Authors

- **Guillermo Aldama-López, MD, PhD** - University Hospital of A Coruña, INIBIC

- **Domingo López-Vázquez, MD** - University Hospital of A Coruña, INIBIC

- **Fernando Rebollal-Leal, MD** - University Hospital of A Coruña, INIBIC

## 🤝 Contributing

We welcome feedback and suggestions. Please:
1. Open an issue for bugs or questions
2. Submit pull requests for improvements
3. Cite our work if you use or extend this analysis

## 📄 License

The content of this repository is dual-licensed:

-   The R code (`analysis.R`) is licensed under the **MIT License**.
-   The data, figures, and all other content are licensed under the **Creative Commons Attribution 4.0 International (CC-BY-4.0) License**.

See the `LICENSE` file for full details.

## 🙏 Acknowledgments

We thank the DIGIT-HF investigators for publishing detailed Kaplan-Meier curves that enabled this secondary analysis. This re-analysis is intended as a constructive methodological contribution to improve the interpretation of clinical trial results.

## 📞 Contact

For questions or collaboration inquiries, please contact:
- **Corresponding author:** Guillermo Aldama-López (guillermo.aldama.lopez@sergas.es)
- **Issues:** Use the [GitHub Issues](https://github.com/galdlop/digit_hf/issues) page

---

**Disclaimer:** This is an independent methodological re-analysis. The authors have no financial or personal conflicts of interest related to the DIGIT-HF trial or digitoxin.

**Last updated:** 18/11/2025

