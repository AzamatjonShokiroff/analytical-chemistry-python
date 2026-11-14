# Analytical Chemistry Data Processing & Method Validation in Python

Python scripts developed for quantitative analytical chemistry workflows, including microplate fluorescence spectroscopy calibration and HPLC raw material standardization.

## Scripts Included
1. **`fluorescence_quenching_calibration.py`**
   - Fits ordinary least squares (OLS) linear regression models to blank-corrected fluorescence quenching data ($\Delta I = I_0 - I$).
   - Computes ICH Q2(R1) Limit of Detection ($\text{LOD} = 3.3 \sigma / S$) and Limit of Quantification ($\text{LOQ} = 10 \sigma / S$).
   - Quantifies unknown samples and exports publication-ready calibration plots.

2. **`hplc_phytochemical_comparison.py`**
   - Evaluates active marker content, precision ($\text{\%RSD}$), and batch-to-batch equivalence across botanical raw material sources using `pandas`, `scipy.stats` (One-Way ANOVA), and `matplotlib`.

## Requirements
```bash
pip install numpy pandas scipy matplotlib
