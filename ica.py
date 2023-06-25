import numpy as np
import pandas as pd
from sklearn.decomposition import FastICA
import matplotlib.pyplot as plt

# Load miRNA expression data for healthy and cancerous samples
healthy_miRNA_data = pd.read_csv("healthy_miRNA_data.csv").values[:, 1:]
cancerous_miRNA_data = pd.read_csv("cancerous_miRNA_data.csv").values[:, 1:]

# Load differentially expressed miRNA data
diff_exp_miRNA_data = pd.read_csv("diff_exp_miRNA_data.csv").values[:, 1:]

# Combine healthy and cancerous miRNA data into one matrix
miRNA_data = np.vstack((healthy_miRNA_data, cancerous_miRNA_data))

# Perform ICA on miRNA data
ica = FastICA(n_components=2)
ica_components = ica.fit_transform(miRNA_data)

# Plot ICA components
plt.scatter(ica_components[:len(healthy_miRNA_data), 0], ica_components[:len(healthy_miRNA_data), 1], c='green', label='Healthy')
plt.scatter(ica_components[len(healthy_miRNA_data):, 0], ica_components[len(healthy_miRNA_data):, 1], c='red', label='Cancerous')
plt.legend()
plt.xlabel("ICA Component 1")
plt.ylabel("ICA Component 2")
plt.show()
