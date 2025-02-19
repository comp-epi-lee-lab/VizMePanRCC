import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path

pickle_path = Path(r'data/pickle_file.pk1')
table = pd.read_pickle(pickle_path)

def compute_significant_cgids(group1, group2, table):
    p_values = {}
    methylation_differences = {}
    
    for column in table.columns[17:]:
        data1 = table.loc[group1, column].dropna()
        data2 = table.loc[group2, column].dropna()
        
        if not data1.empty and not data2.empty:
            _, p_val = stats.mannwhitneyu(data1, data2)
            p_values[column] = p_val
            methylation_differences[column] = abs(data1.mean() - data2.mean())
    
    cgIDs = list(p_values.keys())
    p_vals = list(p_values.values())
    fdr_corrected = stats.false_discovery_control(p_vals)
    
    significant_cgIDs = {cgID for cgID, fdr in zip(cgIDs, fdr_corrected) if (fdr <= 0.05 and methylation_differences[cgID] >= 0.6)}
    
    return significant_cgIDs

AAA = compute_significant_cgids(table["rcc"] == "rcc", table["rcc"] == "normal", table)

BBB = compute_significant_cgids(table["age_at_initial_pathologic_diagnosis"] < 50, 
                                table["age_at_initial_pathologic_diagnosis"] >= 50, table)

XXX = AAA - BBB

young = table["age_at_initial_pathologic_diagnosis"] >= 50
YYY = compute_significant_cgids((table["rcc"] == "rcc") & young, (table["rcc"] == "normal") & young, table)

old = table["age_at_initial_pathologic_diagnosis"] < 50
ZZZ = compute_significant_cgids((table["rcc"] == "rcc") & old, (table["rcc"] == "normal") & old, table)

print("Number of significant cgIDs in AAA:", len(AAA))
print("Number of significant cgIDs in BBB:", len(BBB))
print("Number of significant cgIDs in XXX (AAA - BBB):", len(XXX))
print("Number of significant cgIDs in YYY:", len(YYY))
print("Number of significant cgIDs in ZZZ:", len(ZZZ))