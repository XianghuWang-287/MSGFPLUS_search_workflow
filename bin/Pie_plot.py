import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import matplotlib.pyplot as plt
import numpy as np



cluster_results = pd.read_csv('/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_results/msculster_results/clustering/clusterinfo.tsv', sep='\t')  # Adjust file path and format accordingly


database_results = pd.read_csv('./filtered.tsv', sep='\t')  # Adjust file path and format accordingly
database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
print("Unique values in cluster_results:")
print("Filename:", cluster_results['#Filename'].unique())
print("Scan:", cluster_results['#Scan'].unique())
print("Charge:", cluster_results['#Charge'].unique())


print("\nUnique values in database_results:")
print("MzIDFileName:", database_results['MzIDFileName'].unique())
print("ScanNumber:", database_results['ScanNumber'].unique())
print("Charge:", database_results['Charge'].unique())

unique_spectra_db = database_results[['MzIDFileName', 'ScanNumber']].drop_duplicates()

# Filter cluster_results to only include spectra present in database results
filtered_cluster_results = cluster_results.merge(unique_spectra_db,
                                                  left_on=['#Filename', '#Scan'],
                                                  right_on=['MzIDFileName', 'ScanNumber'],
                                                  how='inner')

# Merge the filtered clusters with the corresponding peptide identifications based on filename, scan number, and charge
merged_data = pd.merge(filtered_cluster_results, database_results,
                       left_on=['#Filename', '#Scan', '#Charge'],
                       right_on=['MzIDFileName', 'ScanNumber', 'Charge'],
                       how='inner')

cluster_peptide_portion = merged_data.groupby('#ClusterIdx')['Peptide'].value_counts(normalize=True)

# Unique peptide types in the dataset
unique_peptide_types = merged_data['Peptide'].unique()

unique_cluster_indices = cluster_peptide_portion.index.get_level_values(0).unique()

# Take the first 100 unique cluster indices
first_100_clusters = unique_cluster_indices[:100]

num_cols = 5

num_clusters = len(first_100_clusters)
num_rows = int(np.ceil(num_clusters / num_cols))

# Create subplots
fig, axes = plt.subplots(num_rows, num_cols, figsize=(30, 5 * num_rows))

# Flatten the axes array and iterate through the first 100 clusters to plot the pie chart
for i, cluster_idx in enumerate(first_100_clusters):
    portion_data = cluster_peptide_portion.loc[cluster_idx]
    unique_peptide_types_cluster = portion_data.index
    row = i // num_cols
    col = i % num_cols
    ax = axes[row, col]
    ax.pie(portion_data, labels=unique_peptide_types_cluster, autopct='%1.1f%%')
    ax.set_title(f'Cluster {cluster_idx} Peptide Fraction')

# Hide empty subplots if any
for i in range(num_clusters, num_rows * num_cols):
    fig.delaxes(axes.flatten()[i])

plt.tight_layout()
plt.show()

