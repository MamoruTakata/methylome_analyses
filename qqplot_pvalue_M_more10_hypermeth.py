import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import probplot

# Load data
df = pd.read_csv('ratio_corr_pvalue_M_more10_hypermeth.bed', sep='\t', header=None)

# Preprocess data
df = df.dropna()  # Remove missing values
p_values = df[4]  # Extract the column with p-values

# Transform p-values to -log10 scale
observed = -np.log10(np.sort(p_values))  # Observed -log10(p-values)
expected = -np.log10(np.linspace(1 / len(p_values), 1, len(p_values)))  # Expected -log10(p-values) under null hypothesis

# Create QQ plot
plt.figure(figsize=(8, 8))
plt.scatter(expected, observed, c="blue", s=10, label="P-values")
plt.plot([0, max(expected)], [0, max(expected)], color="red", linestyle="--", label="Expected")  # Diagonal line (y = x)
plt.xlabel("Expected -log10(P)")
plt.ylabel("Observed -log10(P)")
plt.title("QQ Plot for EWAS")
plt.legend()
plt.tight_layout()

# Save the plot
plt.savefig("qqplot_pvalue_M_4_more10_hypermeth.png")
plt.show()
