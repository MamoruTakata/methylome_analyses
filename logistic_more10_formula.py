import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import glm

input_file = "methylBase_CG_2.txt" 

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

    methyl_cov = opk_methyl_cov + ypk_methyl_cov
    conditions = [0] * len(opk_methyl_cov) + [1] * len(ypk_methyl_cov)

    valid_indices = [i for i, (mC, cov) in enumerate(methyl_cov) if not np.isnan(mC) and not np.isnan(cov)]
    methyl_cov = [methyl_cov[i] for i in valid_indices]
    conditions = [conditions[i] for i in valid_indices]

    opk_count = conditions.count(0)
    ypk_count = conditions.count(1)

    if opk_count > 0 and ypk_count > 0:
        valid_data = pd.DataFrame({
            'mC': [mC for mC, cov in methyl_cov],
            'unmC': [cov - mC for mC, cov in methyl_cov],
            'condition': conditions
        })

        formula = "mC + unmC ~ condition"
        try:
            result = glm(formula=formula, data=valid_data, family=sm.families.Binomial()).fit()
            p_value = result.pvalues["condition"] 
            coefficient = result.params["condition"] 
        except Exception as e:
            p_value = np.nan 
            coefficient = np.nan
    else:
        p_value = np.nan
        coefficient = np.nan  

    output_data.append([chromosome, position, coefficient, p_value])

output_df = pd.DataFrame(output_data, columns=["Chromosome", "Position", "Coefficient", "P_value"])

output_file = "logistic_regression_more10_formula.txt"
output_df.to_csv(output_file, sep="\t", index=False)
print(f"Results with coefficients and filtering saved to {output_file}")
