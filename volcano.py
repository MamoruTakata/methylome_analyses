# Scientific libraries
import numpy as np
import pandas as pd

# Graphic libraries
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set parameters
rcParams['font.size'] = 16
rcParams['figure.figsize'] = (18, 9)

# Generate Figure & Axes instances
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
colors = ["#9eeb30","#b3734d","#888888"]
labels = ["hypermethylated sites in YPK", "hypermethylated sites in OPK", "n.s."]

df = pd.read_table("logistic_regression_more10_formula_qvalue_with_methdiff.txt")

df = df[df["qvalue_slim"] > 0]
down = df[(df["qvalue_slim"] < 10**-2) \
          & (df["Meth_diff"] < -25)]
up = df[(df["qvalue_slim"] < 10**-2) \
        & (df["Meth_diff"] > 25)]
non_sig = df[(df["qvalue_slim"] >= 10**-2) \
             | ((df["Meth_diff"] <= 25) & (df["Meth_diff"] >= -25))]

for i, met in enumerate([down, up, non_sig]):
    ax.scatter(met["Meth_diff"],
               -np.log10(met["qvalue_slim"]),
               label=labels[i],
               color=colors[i],
               alpha=0.5)

# X-axis Settings
# ax.set_xlim(0, 100)
ax.set_xlabel("Difference of methylation level (%)")

# Y-axis Settings
ax.set_ylabel("$-\log_10($q-value$)$")

# Legend Settings

ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

filename = "volcano_compress"
out_dir = "."
for fmt in ["pdf","png", "svg"]:
    plt.savefig("{}/{}.{}".format(out_dir,
                                  filename,
                                  fmt),
                format=fmt,
                bbox_inches="tight")