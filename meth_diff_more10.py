import pandas as pd
import numpy as np

input_file = "methylBase_CG_2.txt" 
meth_diff_output = "meth_diff_more10_results.txt" 

data = pd.read_csv(input_file, sep="\t", comment="#", header=None)

output_data = []

opk_columns = [(i, i + 1) for i in range(4, 37, 3)] 
ypk_columns = [(i, i + 1) for i in range(37, 61, 3)] 

for index, row in data.iterrows():
    chromosome = row[0]
    position = row[1]

    opk_methyl_cov = [
        (row[meth], row[cov]) if row[cov] >= 10 else (np.nan, np.nan)
        for cov, meth in opk_columns
    ]
    ypk_methyl_cov = [
        (row[meth], row[cov]) if row[cov] >= 10 else (np.nan, np.nan)
        for cov, meth in ypk_columns
    ]

    opk_methyl_rates = [mC / cov * 100 for mC, cov in opk_methyl_cov if not np.isnan(mC) and not np.isnan(cov)]
    ypk_methyl_rates = [mC / cov * 100 for mC, cov in ypk_methyl_cov if not np.isnan(mC) and not np.isnan(cov)]

    if len(opk_methyl_rates) > 0 and len(ypk_methyl_rates) > 0:
        mean_opk = np.mean(opk_methyl_rates)
        mean_ypk = np.mean(ypk_methyl_rates)
        meth_diff = mean_opk - mean_ypk
    else:
        meth_diff = np.nan

    output_data.append([chromosome, position, meth_diff])

meth_diff_df = pd.DataFrame(output_data, columns=["Chromosome", "Position", "Meth_diff"])

meth_diff_df.to_csv(meth_diff_output, sep="\t", index=False)
print(f"meth.diff values saved to {meth_diff_output}")