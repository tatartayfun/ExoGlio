import pandas as pd
import numpy as np
import io
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns

differentially_expressed_miRNA = []

with open('GSE112462_series_matrix.txt', 'r') as f:
    data = f.read()

# Find the index of the start and end markers
start_index = data.index("!series_matrix_table_begin")
end_index = data.index("!series_matrix_table_end")

# Extract the data between the start and end markers
extracted_data = data[start_index:end_index]

# Create a DataFrame from the extracted data, setting the first row as the column names
df = pd.read_csv(io.StringIO(extracted_data), sep='\t', header=1, error_bad_lines=False)

# Create a list of the column names to group as "CTRL"
ctrl_columns = ["GSM3070465", "GSM3070466", "GSM3070467", "GSM3070468", "GSM3070469", "GSM3070470", "GSM3070471", "GSM3070472"]

# Create a list of the column names to group as "OLIGO"
oligo_columns = ["GSM3070482", "GSM3070483", "GSM3070484", "GSM3070485", "GSM3070486", "GSM3070487", "GSM3070488", "GSM3070489", "GSM3070490"]

# Select the data for the ctrl_columns and oligo_columns
miR_data = df['ID_REF']
ctrl_data = df[ctrl_columns]
oligo_data = df[oligo_columns]

# Iterate through each miRNA
for i in range(len(miR_data)):
    # Calculate the average expression level in the control group
    ctrl_mean = ctrl_data.iloc[i].mean()
    
    # Calculate the average expression level in the experimental group
    oligo_mean = oligo_data.iloc[i].mean()
    
    # Calculate the fold change
    fold_change = np.log2(oligo_mean / ctrl_mean)
    
    # Perform a t-test to test the hypothesis that the average expression level in the control group is equal to the average expression level in the experimental group
    t, p = ttest_ind(ctrl_data.iloc[i], oligo_data.iloc[i])
    
    # If the p-value is less than 0.05 and the absolute value of the fold change is greater than 2, the miRNA is differentially expressed
    if p < 0.05 and abs(fold_change) > 1:
        differentially_expressed_miRNA.append((miR_data.iloc[i], p, fold_change))

print(differentially_expressed_miRNA)

# Convert the differentially_expressed_miRNA list to a DataFrame
df = pd.DataFrame(differentially_expressed_miRNA, columns=['miRNA', 'p_value', 'fold_change'])

# Save the DataFrame to a .csv file
df.to_csv('oligo_differentially_expressed_miRNA.csv', index=False)

# Concatenate the ctrl_data and oligo_data DataFrames
data = pd.concat([ctrl_data, oligo_data], axis=1)

# Calculate the Pearson correlation matrix
corr = data.corr()

# Visualize the Pearson correlation matrix as a hierarchical heatmap
sns.clustermap(corr, cmap='RdYlGn', linewidths=.5, figsize=(10, 10))
plt.show()