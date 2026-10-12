from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp


# set up project path
PROJECT_DIR = Path.home() / "cfDNA_fragmentomics"

FRAGMENT_DIR = PROJECT_DIR / "result/fragments"
FEATURE_DIR = PROJECT_DIR  / "result/features/fragment_length"

# retrieve all samples' fragment length
samples = {}

for file in FRAGMENT_DIR.glob("*.fragments.tsv"):
    sample_name = file.stem.split(".")[0]
    
    df = pd.read_csv(file, sep="\t")
    lengths = pd.to_numeric(df["length"], errors="coerce").dropna()
    samples[sample_name] = lengths

# descriptive statistics
rows = []

for sample, length in samples.items():
    rows.append({
        "sample": sample,
        "n_fragments": len(length),
        "mean_bp": length.mean(),
        "SD_bp": length.std(),
        "median_bp": length.median(),
        "Q1_bp": length.quantile(0.25),
        "Q3_bp": length.quantile(0.75),
        "IQR_bp": length.quantile(0.75) - length.quantile(0.25),
        "min_bp": length.min(),
        "max_bp": length.max()
    })

pd.DataFrame(rows).to_csv(f"{FEATURE_DIR}/fragment_length_statistics.tsv", sep="\t", index=False)

# histogram of fragment length distribution
sample_order = ["Sample1", "Sample2", "Sample3", "Sample4"]

plt.figure(figsize=(10, 6))
x_left = min(length.min() for length in samples.values())
x_right = max(length.max() for length in samples.values())

for sample in sample_order:
    length = samples[sample]
    plt.hist(
        length,
        label=sample,
        bins=np.arange(x_left-1, x_right+1, 2),
        density=True,
        histtype="step",
        linewidth=0.5
    )

plt.xlabel("Fragment length (bp)", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.title("Fragment Length Distribution", fontsize=14)
plt.xlim(x_left, x_right)
plt.legend(title="Sample")
plt.tight_layout()
plt.savefig(f"{FEATURE_DIR}/fragment_length_histogram.png", dpi=300)
plt.show()

# ANOVA and post hoc test
names = list(samples.keys())
data = list(samples.values())

# test equal variances
test = ""

stat_l, p_l = stats.levene(*data)
if p_l < 0.05:
    # Kruskal-Wallis
    test = "Kruskal-Wallis"
    stat, p = stats.kruskal(*data)
    
    with open(f"{FEATURE_DIR}/kruskal_wallis_result.txt", "w") as f:
        f.write(f"H = {stat:.3f}\n")
        f.write(f"p-value = {p:.4e}\n")
        if p < 0.05:
            # Dunn's test (Holm-adjusted)
            length_arr = np.concatenate(data)
            length_id = np.concatenate([np.repeat(name, len(values)) for name, values in zip(names, data)])
            
            dunn_df = pd.DataFrame({"fragment_length": length_arr, "sample": length_id})
            dunn_result = sp.posthoc_dunn(dunn_df, val_col="fragment_length", group_col="sample", p_adjust="holm")
            
            dunn_result.to_csv(f"{FEATURE_DIR}/dunn_result.tsv", sep="\t")
else:
    # ANOVA
    test = "One-way ANOVA"
    f_stat, p = stats.f_oneway(*data)
    
    with open(f"{FEATURE_DIR}/anova_result.txt", "w") as f:
        f.write(f"f = {f_stat:.3f}")
        f.write(f"p-value = {p:.4e}")
        if p < 0.05:
            # Tukey HSD
            length_arr = np.concatenate(data)
            length_id = np.concatenate([np.repeat(name, len(values)) for name, values in samples.items()])
            
            tukey = pairwise_tukeyhsd(length_arr, length_id, alpha=0.05)
            tukey_result = pd.DataFrame(tukey._results_table.data[1:], columns=tukey._results_table.data[0])
            
            tukey_result.to_csv(f"{FEATURE_DIR}/tukey_result.tsv", sep="\t", index=False)

# box plot
means = [values.mean() for values in data]
plt.figure(figsize=(10, 6))

plt.boxplot(data, tick_labels=names, showfliers=False)
plt.scatter(range(1, len(names) + 1), means, marker="D")

plt.xlabel("Sample")
plt.ylabel("Fragment length (bp)")
plt.title(f"{test}\n(p-value = {p:.4e})")
plt.tight_layout()
plt.savefig(f"{FEATURE_DIR}/fragment_length_boxplot.png", dpi=300)
plt.close()
