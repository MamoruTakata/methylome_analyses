import csv
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

# Set working directory and input file name

file_name = "methylBase_CG_nymph_8.txt"

# Fixed nymph and worker counts for each sample
nymph_worker_counts = {
    'OPK1': (13, 47), 'OPK2': (10, 42), 'OPK6': (3, 18), 'OPK7': (2, 17),
    'OPK8': (6, 29), 'OPK9': (8, 23), 'OPK10': (7, 22), 'OPK11': (6, 24),
    'YPK1': (24, 36), 'YPK2': (11, 20), 'YPK6': (19, 37), 'YPK7': (24, 28),
    'YPK8': (26, 34), 'SK1': (21, 75), 'SK2': (16, 37)
}

# Specify column indices for methylated and unmethylated cytosine counts
methylated_indices = {
    'OPK1': 5, 'OPK2': 8, 'OPK6': 20, 'OPK7': 23, 'OPK8': 26, 'OPK9': 29,
    'OPK10': 32, 'OPK11': 35, 'YPK1': 38, 'YPK2': 41, 'YPK6': 53, 'YPK7': 56,
    'YPK8': 59, 'SK1': 62, 'SK2': 65
}
unmethylated_indices = {
    'OPK1': 6, 'OPK2': 9, 'OPK6': 21, 'OPK7': 24, 'OPK8': 27, 'OPK9': 30,
    'OPK10': 33, 'OPK11': 36, 'YPK1': 39, 'YPK2': 42, 'YPK6': 54, 'YPK7': 57,
    'YPK8': 60, 'SK1': 63, 'SK2': 66
}

# Initialize list to store results
results = []

# Open the methylation data file
with open(file_name, "r") as file:
    f = csv.reader(file, delimiter="\t", skipinitialspace=True)

    for line_number, line in enumerate(f, 1):  # Add line number tracking
        if len(line) < 67:  # Skip lines with insufficient columns
            print(f"Skipping line {line_number}: insufficient columns.")
            continue

        # Get scaffold and locus position
        Scaffold = line[0]
        Locus = line[1]

        # Retrieve methylated and unmethylated cytosine counts for each individual
        try:
            methylated_C = [float(line[methylated_indices[indiv]]) for indiv in methylated_indices]
            unmethylated_C = [float(line[unmethylated_indices[indiv]]) for indiv in unmethylated_indices]
        except IndexError:
            print(f"Error reading methylated/unmethylated columns on line {line_number}.")
            continue

        # Calculate total coverage
        coverages = [m + u for m, u in zip(methylated_C, unmethylated_C)]

        # Exclude data with coverage less than 10
        valid_indices = [i for i, cov in enumerate(coverages) if cov >= 10]
        if not valid_indices:
            print(f"Skipping line {line_number}: all coverages are below 10.")
            continue

        # Extract only valid data
        methylated_C = [methylated_C[i] for i in valid_indices]
        unmethylated_C = [unmethylated_C[i] for i in valid_indices]
        valid_columns = [list(nymph_worker_counts.keys())[i] for i in valid_indices]

        # Calculate M values (log2((m + 1) / (u + 1)))
        try:
            M_values = [np.log2((m + 1) / (u + 1)) for m, u in zip(methylated_C, unmethylated_C)]
        except ZeroDivisionError:
            print(f"Skipping line {line_number}: zero unmethylated counts.")
            continue

        # Skip lines with extreme M values (absolute value > 10)
        if any(np.abs(M) > 10 for M in M_values):
            print(f"Skipping line {line_number} with extreme M-values.")
            continue

        # Calculate nymph ratios (nymph / (nymph + worker))
        nymph_counts = [nymph_worker_counts[col][0] for col in valid_columns]
        worker_counts = [nymph_worker_counts[col][1] for col in valid_columns]
        total_counts = np.array(nymph_counts) + np.array(worker_counts)
        nymph_ratios = np.array(nymph_counts) / total_counts

        try:
            # Calculate Pearson correlation
            correlation, p_value = pearsonr(M_values, nymph_ratios)

            # Append result
            results.append([Scaffold, Locus, correlation, p_value])

        except Exception as e:
            print(f"  Error calculating Pearson correlation on line {line_number}: {e}")
            continue

# Convert results to DataFrame and save as a tab-separated text file
results_df = pd.DataFrame(results, columns=['Scaffold', 'Position', 'Correlation', 'P_value'])
results_df.to_csv("ratio_corr_pvalue_M_more10.out", sep='\t', index=False)

print("Results (Scaffold, Position, Correlation, P_value) have been saved to 'ratio_corr_pvalue_M_more10.out' as a tab-separated file.")