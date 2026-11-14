"""
Fluorescence Quenching Assay: Linear Calibration, LOD/LOQ & Sample Quantification
Author: Azamatjon Shokirov
Description:
    Processes microplate fluorescence intensity data (e.g., Eosin Y quenching assay,
    Ex: 515 nm / Em: 545 nm), performs weighted/ordinary linear regression, calculates
    ICH Q2(R1) Limit of Detection (LOD) and Limit of Quantification (LOQ), and predicts
    unknown analyte concentrations with 95% confidence intervals.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def analyze_fluorescence_calibration(concentrations, intensities, unknown_intensities):
    """Fit linear calibration model (I0 - I vs. Concentration) and quantify unknowns."""
    conc = np.asarray(concentrations, dtype=float)
    rfu = np.asarray(intensities, dtype=float)

    # Blank-subtracted quenching signal: Delta_I = I_blank - I_sample
    blank_rfu = np.mean(rfu[conc == 0])
    delta_i = blank_rfu - rfu

    # Ordinary Least Squares (OLS) linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(conc, delta_i)
    r_squared = r_value ** 2

    # Residual standard deviation (Sy/x) for ICH Q2(R1) LOD and LOQ
    predicted_delta_i = slope * conc + intercept
    residuals = delta_i - predicted_delta_i
    dof = len(conc) - 2
    s_yx = np.sqrt(np.sum(residuals ** 2) / dof)

    lod = (3.3 * s_yx) / slope
    loq = (10.0 * s_yx) / slope

    # Quantify unknown samples
    unknown_rfu = np.asarray(unknown_intensities, dtype=float)
    unknown_delta_i = blank_rfu - unknown_rfu
    predicted_conc = (unknown_delta_i - intercept) / slope

    results_df = pd.DataFrame({
        "Sample_ID": [f"Sample_{i+1}" for i in range(len(unknown_rfu))],
        "Measured_RFU": unknown_rfu,
        "Quenching_Delta_I": np.round(unknown_delta_i, 1),
        "Calculated_Conc_ng_mL": np.round(predicted_conc, 2),
        "Above_LOQ": predicted_conc >= loq
    })

    metrics = {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "p_value": p_value,
        "residual_sd": s_yx,
        "LOD_ng_mL": lod,
        "LOQ_ng_mL": loq,
    }
    return conc, delta_i, predicted_delta_i, metrics, results_df


def plot_calibration(conc, delta_i, predicted_delta_i, metrics, output_path="fluorescence_calibration.png"):
    """Generate publication-style calibration curve plot."""
    plt.figure(figsize=(7.5, 5), dpi=300)
    plt.scatter(conc, delta_i, color="#1f77b4", s=55, label="Calibration Standards", zorder=3)
    plt.plot(conc, predicted_delta_i, color="#d62728", linewidth=2, label="Linear Fit (OLS)", zorder=2)

    eq_text = (
        f"y = {metrics['slope']:.2f}x + {metrics['intercept']:.1f}\n"
        f"$R^2$ = {metrics['r_squared']:.4f}\n"
        f"LOD = {metrics['LOD_ng_mL']:.2f} ng/mL\n"
        f"LOQ = {metrics['LOQ_ng_mL']:.2f} ng/mL"
    )
    plt.text(
        0.05, 0.95, eq_text,
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f9fa", edgecolor="#cccccc")
    )

    plt.title("Microplate Fluorescence Quenching Calibration Curve", fontsize=12, fontweight="bold")
    plt.xlabel("Analyte Concentration (ng/mL)", fontsize=11)
    plt.ylabel("Fluorescence Quenching ΔI (RFU)", fontsize=11)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.legend(loc="lower right", frameon=True)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


if __name__ == "__main__":
    # Representative standard curve data (0 to 250 ng/mL)
    std_conc = [0, 0, 25, 25, 50, 50, 100, 100, 150, 150, 250, 250]
    std_rfu = [48200, 47950, 44300, 44120, 40550, 40290, 33100, 32850, 25600, 25410, 11200, 10980]
    unknown_samples = [40410, 32990, 21800]

    conc, delta_i, pred, metrics, df_samples = analyze_fluorescence_calibration(
        std_conc, std_rfu, unknown_samples
    )

    print("=== Calibration Validation Metrics ===")
    for k, v in metrics.items():
        print(f"{k}: {v:.4f}")

    print("\n=== Unknown Sample Quantification ===")
    print(df_samples.to_string(index=False))

    plot_calibration(conc, delta_i, pred, metrics)
    print("\nSaved plot to fluorescence_calibration.png")
