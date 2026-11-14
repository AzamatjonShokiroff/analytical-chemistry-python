"""
HPLC Phytochemical Standardization & Batch Equivalence Analysis
Author: Azamatjon Shokirov
Description:
    Evaluates batch-to-batch consistency (Mean, SD, %RSD) of active marker compounds
    quantified by HPLC-UV across locally cultivated botanical raw materials versus
    imported commercial reference standards using One-Way ANOVA and Tukey-style comparisons.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def summarize_batches(df):
    """Compute summary statistics and %RSD for each batch group."""
    summary = df.groupby("Source")["Active_Content_mg_g"].agg(
        Replicates="count",
        Mean_mg_g="mean",
        SD_mg_g="std"
    ).reset_index()

    summary["RSD_percent"] = (summary["SD_mg_g"] / summary["Mean_mg_g"]) * 100
    return summary.round(3)


def run_anova(df):
    """Perform One-Way ANOVA across sample sources."""
    groups = [group["Active_Content_mg_g"].values for _, group in df.groupby("Source")]
    f_stat, p_val = stats.f_oneway(*groups)
    return f_stat, p_val


def plot_batch_comparison(df, output_path="hplc_batch_comparison.png"):
    """Generate boxplot + individual replicate overlay for HPLC marker content."""
    sources = df["Source"].unique()
    data = [df[df["Source"] == s]["Active_Content_mg_g"].values for s in sources]

    plt.figure(figsize=(7.5, 5), dpi=300)
    box = plt.boxplot(data, labels=sources, patch_artist=True, widths=0.45)

    colors = ["#2ca02c", "#1f77b4", "#ff7f0e"]
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.35)

    # Overlay individual replicate points
    for idx, values in enumerate(data, start=1):
        jitter = np.random.normal(0, 0.04, size=len(values))
        plt.scatter(np.repeat(idx, len(values)) + jitter, values, color="#222222", s=35, zorder=3)

    plt.title("HPLC Active Marker Quantification Across Raw Material Sources", fontsize=11, fontweight="bold")
    plt.ylabel("Active Marker Content (mg/g dry weight)", fontsize=10)
    plt.grid(True, axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


if __name__ == "__main__":
    np.random.seed(42)
    dataset = pd.DataFrame({
        "Source": (
            ["Local Cultivated (Batch A)"] * 6 +
            ["Local Cultivated (Batch B)"] * 6 +
            ["Imported Commercial Ref"] * 6
        ),
        "Active_Content_mg_g": [
            14.8, 15.1, 14.9, 15.2, 14.7, 15.0,
            14.6, 14.9, 15.0, 14.8, 15.1, 14.7,
            15.0, 15.3, 14.9, 15.2, 15.1, 14.8,
        ]
    })

    summary_table = summarize_batches(dataset)
    f_stat, p_val = run_anova(dataset)

    print("=== HPLC Standardization Summary ===")
    print(summary_table.to_string(index=False))
    print(f"\nOne-Way ANOVA: F = {f_stat:.4f}, p = {p_val:.4f}")

    plot_batch_comparison(dataset)
    print("Saved comparison chart to hplc_batch_comparison.png")
