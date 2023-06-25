#token: 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

import pandas as pd
import numpy as np
import io
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

differentially_expressed_miRNA = []
top_peaks = []
    
# Read the .txt file into a DataFrame
df = pd.read_csv('GSE122488_normalized_microRNA_counts.txt', delimiter='\t', header=3)[2:]

# Create a list of the column names to group as "CTRL"
ctrl_columns = ['GBM1', 'GBM2',
       'GBM3', 'GBM4', 'GBM5', 'GBM6', 'GBM7', 'GBM8', 'GBM9', 'GBM10',
       'GBM11', 'GBM12', 'GII-III.1', 'GII-III.2', 'GII-III.3', 'GII-III.4',
       'GII-III.5', 'GII-III.6', 'GII-III.7', 'GII-III.8', 'GII-III.9',
       'GII-III.10']

# Create a list of the column names to group as "GLIO"
glio_columns = ['HC1', 'HC2', 'HC3', 'HC4', 'HC5', 'HC6', 'HC7', 'HC8', 'HC9',
       'HC10', 'HC11', 'HC12', 'HC13', 'HC14', 'HC15', 'HC16']

# Select the data for the ctrl_columns and glio_columns
miR_data = df['miRNA']
ctrl_data = df[ctrl_columns]
glio_data = df[glio_columns]

# Iterate through each miRNA
for i in range(len(miR_data)):
    # Calculate the average expression level in the control group
    ctrl_mean = ctrl_data.iloc[i].mean()
    
    # Calculate the average expression level in the experimental group
    glio_mean = glio_data.iloc[i].mean()
    
    # Calculate the fold change
    fold_change = np.log2(glio_mean / ctrl_mean)
    
    # Perform a t-test to test the hypothesis that the average expression level in the control group is equal to the average expression level in the experimental group
    t, p = ttest_ind(ctrl_data.iloc[i], glio_data.iloc[i])
    
    # If the p-value is less than 0.05 and the absolute value of the fold change is greater than 2, the miRNA is differentially expressed
    if p < 0.05 and abs(fold_change) > 1:
        differentially_expressed_miRNA.append((miR_data.iloc[i], p, fold_change))
        top_peaks.append(i)  # Store the index of the differentially expressed miRNA

print(len(differentially_expressed_miRNA))

# Convert the differentially_expressed_miRNA list to a DataFrame
df = pd.DataFrame(differentially_expressed_miRNA, columns=['miRNA', 'p_value', 'fold_change'])

# Save the DataFrame to a .csv file
df.to_csv('differentially_expressed_miRNA.csv', index=False)

# Concatenate the ctrl_data and glio_data DataFrames
data = pd.concat([ctrl_data, glio_data], axis=1)

# Calculate the Pearson correlation matrix
corr = data.corr()

# Visualize the Pearson correlation matrix as a hierarchical heatmap
sns.clustermap(corr, cmap='RdYlGn', linewidths=.5, figsize=(10, 10))
plt.show()

# create a DataFrame to store the mz_axis data
top_peaks_df = pd.DataFrame(top_peaks, columns=['peak_idx'])

# save the top_feat_idx_df DataFrame to an Excel file
with pd.ExcelWriter('top_peaks.xlsx') as writer:
    top_peaks_df.to_excel(writer, index=False)
