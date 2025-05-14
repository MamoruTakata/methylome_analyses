import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

file_path = "logistic_regression_more10_formula_qvalue_with_methdiff.txt"  
df = pd.read_csv(file_path, sep="\t")

df["log_qvalue"] = -np.log10(df["qvalue_slim"])

chromosomes = sorted(df["Chromosome"].unique()) 
chrom_colors = {chrom: "black" if i % 2 == 0 else "gray" for i, chrom in enumerate(chromosomes)}

offsets = {chrom: i * 5e7 for i, chrom in enumerate(chromosomes)} 

df["Plot_Position"] = df["Position"] + df["Chromosome"].map(offsets)

opk_high_meth = (df["Meth_diff"] > 25) & (df["qvalue_slim"] <= 0.01)
ypk_high_meth = (df["Meth_diff"] < -25) & (df["qvalue_slim"] <= 0.01)

fig, ax = plt.subplots(figsize=(12, 6))
for chrom in chromosomes:
    subset = df[(df["Chromosome"] == chrom) & (df["Meth_diff"] > 0)]
    ax.scatter(
        subset["Plot_Position"],
        subset["log_qvalue"],
        color=["#b3734d" if opk_high_meth.loc[idx] else chrom_colors[chrom] for idx in subset.index],
        s=10,
    )

ax.set_xlabel("Chromosome", fontsize=12)
ax.set_ylabel("-log10(qvalue)", fontsize=12)
ax.set_title("Manhattan Plot of Methylation Differences (OPK)", fontsize=14)
ax.axhline(0, color="black", linestyle="--", linewidth=1)

ax.set_xticks([offsets[chrom] + 2.5e7 for chrom in chromosomes])
ax.set_xticklabels(chromosomes)

ax.legend().remove()

opk_output_path = "manhattan_plot_OPK.png"
plt.savefig(opk_output_path, dpi=300, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(12, 6))
for chrom in chromosomes:
    subset = df[(df["Chromosome"] == chrom) & (df["Meth_diff"] < 0)]
    ax.scatter(
        subset["Plot_Position"],
        subset["log_qvalue"],
        color=["#9eeb30" if ypk_high_meth.loc[idx] else chrom_colors[chrom] for idx in subset.index],
        s=10,
    )

ax.set_xlabel("Chromosome", fontsize=12)
ax.set_ylabel("-log10(qvalue)", fontsize=12)
ax.set_title("Manhattan Plot of Methylation Differences (YPK)", fontsize=14)
ax.axhline(0, color="black", linestyle="--", linewidth=1)

ax.set_xticks([offsets[chrom] + 2.5e7 for chrom in chromosomes])
ax.set_xticklabels(chromosomes)

ax.legend().remove()

ypk_output_path = "manhattan_plot_YPK.png"
plt.savefig(ypk_output_path, dpi=300, bbox_inches="tight")
plt.close()

opk_output_path, ypk_output_path