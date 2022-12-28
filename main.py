import pandas as pd
import numpy as np
import io
from scipy.stats import f_oneway

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
for i in range(len(df)):
    # Perform one-way ANOVA test
    stat, p = f_oneway(ctrl_data.iloc[i], glio_data.iloc[i])
    
    # If the p-value is less than 0.05, the miRNA is differentially expressed
    if p < 0.05:
        differentially_expressed_miRNA.append((df['ID_REF'].iloc[i], p))

print(differentially_expressed_miRNA)