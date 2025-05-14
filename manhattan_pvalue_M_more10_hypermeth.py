import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv('ratio_corr_pvalue_M_more10_hypermeth.bed', sep='\t', header=None)
df = df.dropna()

# Calculate -log10(P) for P-values
df['pvalue'] = -np.log10(df[4])
df['ind'] = range(len(df))

# P-value threshold
p_threshold = 1e-5
log_p_threshold = -np.log10(p_threshold)

# Initialize colors based on threshold
df['color'] = 'grey'
df['highlight'] = np.where(df['pvalue'] < log_p_threshold, False, True)

# Alternate colors (grey and black) for each chromosome group
colors_chromosomes = ['grey', 'black']
df_grouped = df.groupby(df[0])

for num, (name, group) in enumerate(df_grouped):
    mask = ~group['highlight']  # Apply only to sites above the threshold
    df.loc[group.index[mask], 'color'] = colors_chromosomes[num % 2]

# Prepare plot
fig = plt.figure(figsize=(14, 8))
ax = fig.add_subplot(111)

x_labels = []
x_labels_pos = []

# Plot for each chromosome group
for num, (name, group) in enumerate(df_grouped):
    group.plot(
        kind='scatter',
        x='ind',
        y='pvalue',
        color=group['color'],
        ax=ax,
        s=10,  # Point size
        alpha=0.7  # Point transparency
    )
    # Highlight sites with P < 1e-5 using sky blue
    group[group['highlight']].plot(
        kind='scatter',
        x='ind',
        y='pvalue',
        color='skyblue',
        ax=ax,
        s=10,
        alpha=0.7
    )
    x_labels.append(name)
    x_labels_pos.append(
        (group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2)
    )

# X-axis labels
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)

# Add threshold lines
ax.axhline(y=log_p_threshold, color='skyblue', linestyle='--', linewidth=1.5)
ax.axhline(y=-np.log10(0.05/1631109), color='red', linestyle='--', linewidth=1.5)

# Set axis limits
ax.set_xlim([0, len(df)])
ax.set_ylim([0, 12])

# Axis labels
ax.set_xlabel('Chr')
ax.set_ylabel('-log10(P)')

# Save figure
fig.savefig('manhattan_pvalue_M_4_more10_hypermeth.png', dpi=300)

print("Manhattan plot with chromosome-based coloring and highlights saved to 'manhattan_pvalue_M_4_more10_hypermeth.png'.")