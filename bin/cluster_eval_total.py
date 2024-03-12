import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import matplotlib.pyplot as plt

def calculate_cluster_purity(cluster):
    total_spectra = len(cluster)
    if total_spectra == 1:
        return 1

    # Filter out 'unknown' peptides and calculate the count of identified peptides per file in the cluster
    identified_peptides = cluster[cluster['Peptide'] != 'unknown']
    identified_peptides_counts = identified_peptides.groupby('Peptide').size()
    max_count = identified_peptides_counts.max() if not identified_peptides_counts.empty else 0

    max_fraction = max_count/total_spectra

    return max_fraction

# Load cluster results
cluster_results = pd.read_csv('/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_results/msculster_results/clustering/mscluster_clusterinfo.tsv', sep='\t')  # Adjust file path and format accordingly

# Load database search results
database_results = pd.read_csv('./filtered.tsv', sep='\t')  # Adjust file path and format accordingly
database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
print("Unique values in cluster_results:")
print("Filename:", cluster_results['#Filename'].unique())
print("ClusterIdx:", cluster_results['#ClusterIdx'].unique())


# Create a merged dataset with left join to include all cluster results and match with database results
merged_data = cluster_results.merge(database_results, left_on=['#Filename', '#Scan'], right_on=['MzIDFileName', 'ScanNumber'], how='left')

# Fill 'Peptide' column with 'unknown' for unidentified spectra
merged_data['Peptide'].fillna('unknown', inplace=True)
print("Unique values in merged_data:")
print("Filename:", merged_data['#Filename'].unique())
print("ClusterIdx:", merged_data['#ClusterIdx'].unique())
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