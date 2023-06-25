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
top_feat_idx = top_feat_idx[:10]
    
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
