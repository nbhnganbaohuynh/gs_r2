from pathlib import Path
from itertools import product
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# set up project path
PROJECT_DIR = Path.home() / "cfDNA_fragmentomics"

FRAGMENT_DIR = PROJECT_DIR / "result/fragments"
FEATURE_DIR = PROJECT_DIR  / "result/features/end_motif"

# all possible 4-mer motifs
MOTIFS = ["".join(x) for x in product("ACGT", repeat=4)]

# motif count, frequency, and MDS
motif_count = pd.DataFrame(index=MOTIFS)
motif_count.index.name = "end_motif"

motif_freq = pd.DataFrame(index=MOTIFS)
motif_freq.index.name = "end_motif"

sample_order = ["Sample1", "Sample2", "Sample3", "Sample4"]
sample_mds = pd.DataFrame(index=sample_order)
sample_mds.index.name = "sample"

for file in FRAGMENT_DIR.glob("*.fragments.tsv"):
    sample_name = file.stem.split(".")[0]
    
    df = pd.read_csv(file, sep="\t")
    end_motifs = pd.concat([df["left_motif"], df["right_motif"]], ignore_index=True).astype("string").str.upper()
    end_motifs = end_motifs[end_motifs.str.fullmatch(r"[ACGT]{4}", na=False)]
    
    # count occurrences of motifs
    count = end_motifs.value_counts().reindex(MOTIFS, fill_value=0).astype(int)
    total = count.sum()
    freq = count / total
    
    motif_count[sample_name] = count
    motif_freq[sample_name] = freq
    
    # calculate MDS
    p = freq[freq > 0].to_numpy()
    mds = (-np.sum(p * np.log(p)) / np.log(256))
    
    sample_mds.loc[sample_name, "MDS"] = mds

# save frequency and MDS
motif_freq.to_csv(f"{FEATURE_DIR}/end_motif_frequencies.tsv", sep="\t")
sample_mds.to_csv(f"{FEATURE_DIR}/end_motif_mds.tsv", sep="\t")

# top 20 dominant motifs of each sample
fig, axes = plt.subplots(len(sample_order), 1, figsize=(12, 3.5 * len(sample_order)), squeeze=False)

for ax, sample in zip(axes[:, 0], sample_order):
    top = motif_freq[sample].nlargest(20)

    ax.bar(top.index, top.values * 100)
    ax.set_title(sample)
    ax.set_ylabel("Frequency (%)")
    ax.tick_params(axis="x", rotation=60)

axes[-1, 0].set_xlabel("End motif")
fig.tight_layout()
fig.savefig(f"{FEATURE_DIR}/end_motif_top20.png", dpi=300)
plt.close(fig)

# motif starts with "CC"
cc_motifs = [motif for motif in MOTIFS if motif.startswith("CC")]
plot_df = motif_freq.loc[cc_motifs, sample_order] * 100

x = np.arange(len(cc_motifs))
width = 0.2

plt.figure(figsize=(14, 6))

for i, sample in enumerate(sample_order):
    plt.bar(x + (i - 1.5) * width, plot_df[sample], width=width, label=sample)

plt.xticks(x, cc_motifs, rotation=45)
plt.xlabel("End motif")
plt.ylabel("Frequency")
plt.title('Frequencies of End Motifs Starting with "CC"')
plt.legend(title="Sample")
plt.grid(axis="y", linestyle="--", alpha=0.3)
plt.tight_layout()
plt.savefig(f"{FEATURE_DIR}/cc_motif_frequencies.png", dpi=300, bbox_inches="tight")
plt.close()

# dot plot of MDS
plt.figure(figsize=(10, 6))

plt.scatter(sample_mds.index, sample_mds["MDS"], s=100)

for sample, row in sample_mds.iterrows():
    plt.annotate(
        f'{row["MDS"]:.4f}',
        (sample, row["MDS"]),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center"
    )

plt.ylabel("MDS")
plt.ylim(
    max(0, sample_mds["MDS"].min() - 0.1),
    min(1.0, sample_mds["MDS"].max() + 0.1)
)
plt.title("End-motif diversity by sample")
plt.tight_layout()
plt.savefig(f"{FEATURE_DIR}/end_motif_mds.png", dpi=300)
plt.close()

# compare frequency between each two samples pair
freq_diff = pd.DataFrame(columns=[
    "sample_1",
    "sample_2",
    "motif",
    "frequency_1",
    "frequency_2",
    "difference_in_percentages",
])

for i in range(len(sample_order)):
    for j in range(i + 1, len(sample_order)):
        s1, s2 = sample_order[i], sample_order[j]
        
        for motif in MOTIFS:
            f1 = motif_freq.loc[motif, s1]
            f2 = motif_freq.loc[motif, s2]
            
            diff = (f1 - f2) * 100
            freq_diff.loc[len(freq_diff)] = [s1, s2, motif, f1, f2, diff]

freq_diff["abs_difference"] = (freq_diff["difference_in_percentages"].abs())
freq_diff = freq_diff.sort_values(["sample_1", "sample_2", "abs_difference"], ascending=[True, True, False])

freq_diff.to_csv(f"{FEATURE_DIR}/end_motif_differential.tsv", sep="\t", index=False)

# plot top 5 end motifs with biggest differences between each samples pair
top5_df = (
    freq_diff.sort_values("abs_difference", ascending=False)
    .groupby(["sample_1", "sample_2"], sort=True)
    .head(5)
    .copy()
)

pairs = list(top5_df.groupby(["sample_1", "sample_2"], sort=True))

fig, axes = plt.subplots(3, 2, figsize=(15, 16), constrained_layout=True)
axes = axes.flatten()

for ax, ((s1, s2), group) in zip(axes, pairs):
    group = group.sort_values("difference_in_percentages")
    
    ax.barh(group["motif"], group["difference_in_percentages"])
    ax.axvline(0, linewidth=1, color="black")
    ax.set_title(f"{s1} vs {s2}")
    ax.set_xlabel("Frequency difference (%)")
    ax.set_ylabel("End motif")

fig.suptitle("Top 5 End-Motif Frequency Differences per Sample Pair", fontsize=15)

fig.savefig(f"{FEATURE_DIR}/top5_freq_diff_all_pairs.png", dpi=300, bbox_inches="tight")
plt.close(fig)
