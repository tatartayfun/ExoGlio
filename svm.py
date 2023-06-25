import pandas as pd
import numpy as np
import io
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA, FastICA
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, recall_score, precision_score, accuracy_score, confusion_matrix
import seaborn as sns
from skopt import BayesSearchCV
from sklearn.model_selection import GridSearchCV

differentially_expressed_miRNA = []

# specify the path to the Excel file containing mz_axis
peaks_file = 'top_peaks.xlsx'

# read the mz_axis file and store it in a numpy array
top_peaks_df = pd.read_excel(peaks_file, header=0)
top_feat_idx = top_peaks_df.to_numpy()
top_feat_idx = top_feat_idx[:25]
print(top_feat_idx[22:23])
#top_feat_idx = top_feat_idx[22:23]

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

# select only specific rows based on top_feat_idx
miR_data_sel = miR_data.iloc[np.transpose(top_feat_idx)[0]]
ctrl_data_sel = ctrl_data.iloc[np.transpose(top_feat_idx)[0]]
glio_data_sel = glio_data.iloc[np.transpose(top_feat_idx)[0]]

# Concatenate the dataframes along the column axis (axis=1)
combined_data = pd.concat([ctrl_data_sel, glio_data_sel], axis=1).to_numpy().T
# merge the selected dataframes and create a label array
y = np.hstack([np.ones(len(ctrl_data_sel.columns)), np.full(len(glio_data_sel.columns), 2)])

X_train, X_test, y_train, y_test = train_test_split(combined_data, y, test_size=0.3, random_state=42, stratify=y)

svm = SVC(kernel='linear', C = 5.44, gamma=1.0)
"""
# define the parameter search space for the SVM classifier
param_space = {
    'C': (0.01, 10.0, 'log-uniform'),
    'kernel': ['linear', 'rbf'],
    'gamma': (0.001, 1.0, 'log-uniform')
}

# perform Bayesian optimization
opt = BayesSearchCV(
    SVC(),
    param_space,
    n_iter=50,
    scoring='precision',
    n_jobs=-1,
    cv=5
)
opt.fit(X_train, y_train)

# Print the best parameters and score found during optimization
print("Best parameters found:")
print(opt.best_params_)
print("Best accuracy found:")
print(opt.best_score_)

# Get the optimized classifier
svm = opt.best_estimator_
"""
svm.fit(X_train, y_train)

y_pred = svm.predict(X_test)

# make predictions on the test set and print the accuracy
recall = recall_score(y_test, y_pred, average='macro')
precision = precision_score(y_test, y_pred, average='macro')
accuracy = accuracy_score(y_test, y_pred)
print(f'Test recall rate: {recall}')
print(f'Test precision: {precision}')
print(f'Test accuracy: {accuracy}')#print(classification_report(y_test, y_pred))

# Visualize the confusion matrix using a heatmap
cm = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', annot_kws={"fontsize": 14})
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix')

plt.show()
