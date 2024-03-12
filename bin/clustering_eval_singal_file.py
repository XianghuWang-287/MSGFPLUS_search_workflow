import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import matplotlib.pyplot as plt

def calculate_cluster_purity(cluster):
    total_spectra = len(cluster)
    if total_spectra == 1:
        return 1

    # Calculate the count of peptides per file in the cluster
    file_peptide_counts = cluster.groupby('#Filename')['Peptide'].apply(lambda x: x.value_counts().max() / len(x))

    max_fraction = file_peptide_counts.max()

    return max_fraction


cluster_results = pd.read_csv('/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_results/msculster_results/clustering/mscluster_clusterinfo.tsv', sep='\t')  # Adjust file path and format accordingly


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


# Calculate the purity for each cluster
cluster_purity = merged_data.groupby('#ClusterIdx').apply(calculate_cluster_purity)

# Calculate the number of spectra in each cluster
cluster_size = merged_data.groupby('#ClusterIdx').size()

purity_by_cluster_size = pd.DataFrame({'ClusterSize': cluster_size, 'ClusterPurity': cluster_purity})

# Group by cluster size and calculate average purity for each group
average_purity_by_cluster_size = purity_by_cluster_size.groupby('ClusterSize')['ClusterPurity'].mean()

# Plot average purity grouped by cluster size
plt.figure(figsize=(10, 6))
plt.plot(average_purity_by_cluster_size.index, average_purity_by_cluster_size.values, marker='o')
# plt.scatter(cluster_size, cluster_purity, alpha=0.5)
plt.xlabel('Cluster Size (Number of Spectra)')
plt.ylabel('Average Cluster Purity')
plt.title('Average Cluster Purity by Cluster Size')
plt.grid(True)
plt.show()