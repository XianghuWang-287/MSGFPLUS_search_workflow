import pandas as pd
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import matplotlib.pyplot as plt
import numpy as np
import copy

def calculate_cluster_purity(cluster):
    total_spectra = len(cluster)
    if total_spectra == 1:
        return 1

    identified_peptides_counts = cluster.groupby('Peptide').size()
    max_count = identified_peptides_counts.max() if not identified_peptides_counts.empty else 0

    max_fraction = max_count/total_spectra

    return max_fraction

def mscluster_purity(cluster_results,database_results):
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzML')
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
    # Create a merged dataset with left join to include all cluster results and match with database results
    merged_data = cluster_results.merge(database_results, left_on=['#Filename', '#Scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    merged_data = merged_data[merged_data['Peptide'] != 'unknown']
    return merged_data.groupby('#ClusterIdx').apply(calculate_cluster_purity),merged_data.groupby('#ClusterIdx').size()

def falcon_purity(cluster_results,database_results):
    # Create a merged dataset with left join to include all cluster results and match with database results
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    largest_cluster_index = merged_data['cluster'].value_counts().idxmax()

    # Print the cluster index of the largest cluster
    print("Cluster Index of the Largest Cluster:", largest_cluster_index)
    return merged_data.groupby('cluster').apply(calculate_cluster_purity), merged_data.groupby('cluster').size()

def maracluster_purity(cluster_results,database_results):
    database_results['MzIDFileName'] = './mzML/' + database_results['MzIDFileName']
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')
    return merged_data.groupby('cluster').apply(calculate_cluster_purity), merged_data.groupby('cluster').size()

def assign_size_range_bin(size):
    if size == 1:
        return '1'
    else:
        upper_limit = 2
        while size > upper_limit:
            upper_limit *= 2
        return f"{upper_limit//2+1} to {upper_limit}"



# Load cluster results
mscluster_results = pd.read_csv('../data/Combine_results/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
falcon_results = pd.read_csv('../data/Combine_results/Falcon_cluster_info.tsv',sep='\t')  # Adjust file path and format accordingly
database_results = pd.read_csv('./Combine_test_filtered.tsv', sep='\t')  # Adjust file path and format accordingly
maracluster_results= pd.read_csv('../data/Combine_results/MaRaCluster_processed.clusters_p10.tsv', sep='\t')




mscluster_purity,mscluster_size = mscluster_purity(mscluster_results,copy.copy(database_results))
falcon_purity,falcon_size =falcon_purity(falcon_results,copy.copy(database_results))
maracluster_purity,maracluster_size = maracluster_purity(maracluster_results,copy.copy(database_results))

purity_by_cluster_size_mscluster = pd.DataFrame({'ClusterSize': mscluster_size, 'ClusterPurity': mscluster_purity})
purity_by_cluster_size_mscluster['Bin'] = purity_by_cluster_size_mscluster['ClusterSize'].apply(assign_size_range_bin)
# Group by cluster size and calculate average purity for each group
average_purity_by_bin_mscluster = purity_by_cluster_size_mscluster.groupby('Bin')['ClusterPurity'].mean().reset_index()

# Since the bins are non-numeric and custom, we'll sort them by the lower bound of each range for plotting
average_purity_by_bin_mscluster['SortKey'] = average_purity_by_bin_mscluster['Bin'].apply(lambda x: int(x.split(' ')[0]))
average_purity_by_bin_mscluster = average_purity_by_bin_mscluster.sort_values('SortKey')

purity_by_cluster_size_falcon = pd.DataFrame({'ClusterSize':falcon_size, 'ClusterPurity':falcon_purity})
purity_by_cluster_size_falcon['Bin'] = purity_by_cluster_size_falcon['ClusterSize'].apply(assign_size_range_bin)
average_purity_by_bin_falcon = purity_by_cluster_size_falcon.groupby('Bin')['ClusterPurity'].mean().reset_index()
average_purity_by_bin_falcon['SortKey'] = average_purity_by_bin_falcon['Bin'].apply(lambda x: int(x.split(' ')[0]))
average_purity_by_bin_falcon = average_purity_by_bin_falcon.sort_values('SortKey')

purity_by_cluster_size_maracluster = pd.DataFrame(
    {'ClusterSize': maracluster_size, 'ClusterPurity': maracluster_purity})
purity_by_cluster_size_maracluster['Bin'] = purity_by_cluster_size_maracluster['ClusterSize'].apply(
    assign_size_range_bin)
average_purity_by_bin_maracluster = purity_by_cluster_size_maracluster.groupby('Bin')[
    'ClusterPurity'].mean().reset_index()
average_purity_by_bin_maracluster['SortKey'] = average_purity_by_bin_maracluster['Bin'].apply(
    lambda x: int(x.split(' ')[0]))
average_purity_by_bin_maracluster = average_purity_by_bin_maracluster.sort_values('SortKey')

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(average_purity_by_bin_mscluster['Bin'], average_purity_by_bin_mscluster['ClusterPurity'], marker='o', linestyle='-',label='mscluster')
plt.plot(average_purity_by_bin_falcon['Bin'],average_purity_by_bin_falcon['ClusterPurity'], marker='o', linestyle='-',label='falcon')
plt.plot(average_purity_by_bin_maracluster['Bin'],
         average_purity_by_bin_maracluster['ClusterPurity'],
         marker='o', linestyle='-', label='maracluster')
plt.xlabel('Cluster Size Range')
plt.ylabel('Average Purity')
plt.title('Average Cluster Purity by Size Range')
plt.ylim([0,1.1])
plt.xticks(rotation=45)
plt.grid(True, which="both", ls="--", linewidth=0.5)
plt.legend()
plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels
plt.show()