import pandas as pd
import numpy as np
import io
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

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

# Create a list of the column names to group as "GLIO"
glio_columns = ["GSM3070491", "GSM3070492", "GSM3070493", "GSM3070494", "GSM3070495", "GSM3070496", "GSM3070497", "GSM3070498", "GSM3070499", "GSM3070500"]

# Select the data for the ctrl_columns and glio_columns
miR_data = df['ID_REF']
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

print(differentially_expressed_miRNA)

# Select the differentially expressed miRNA data
differentially_expressed_miRNA_data = df.loc[:, ctrl_columns + glio_columns]

# Perform PCA
pca = PCA(n_components=2)
pca.fit(differentially_expressed_miRNA_data)

# Get the transformed data
transformed_data = pca.transform(differentially_expressed_miRNA_data)

# Get the first and second principal components
pc1 = transformed_data[:, 0]
pc2 = transformed_data[:, 1]

# Iterate through each row of the data
for i in range(len(pc1)):
    # If the data is from the CTRL group, plot it with a blue dot
    if i < len(ctrl_columns):
        plt.scatter(pc1[i], pc2[i], color='blue')
    # If the data is from the GLIO group, plot it with a red dot
    else:
        plt.scatter(pc1[i], pc2[i], color='red')

# Add labels and a title to the plot
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.title('PCA of Differentially Expressed miRNA Data')

# Show the plot
plt.show()