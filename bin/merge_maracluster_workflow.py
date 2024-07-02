import pandas as pd
import os
import sys
import argparse
import pymzml
def process_maracluster_results(input_file):
    # Read the input file, considering the format described
    # Each spectrum is on a new line, clusters are separated by an empty line
    # Columns are separated by tabs
    # Use on_bad_lines='skip' to skip any malformed lines
    df = pd.read_csv(input_file, sep='\t', header=None, names=['filename', 'scan', 'cluster'], on_bad_lines='skip')

    # Write the processed DataFrame to a new TSV file
    #df.to_csv(output_file, sep='\t', index=False)
    return df


def enrich_cluster_results(cluster_results,mzml_folder_path):
    cluster_results['filename'] = cluster_results['filename'].apply(os.path.basename).apply(lambda x: os.path.splitext(x)[0])
    print(cluster_results.head())
    workflow_results = pd.read_csv(mzml_folder_path)
    workflow_results['original_path'] = workflow_results['original_path'].apply(os.path.basename).apply(lambda x: os.path.splitext(x)[0])
    print(workflow_results.head())
    cluster_results = cluster_results.merge(workflow_results, left_on=['filename', 'scan'],
                                        right_on=['original_path', 'scan'], how='inner')
    print(max(cluster_results ['retention_time']))

    return cluster_results

# Specify the input and output file paths
input_file = '/home/user/LabData/XianghuData/MS_Cluster_datasets/PXD021518_convert/maracluster_output/MaRaCluster.clusters_p5.tsv'
input_workflow_reuslt='/home/user/LabData/XianghuData/MS_Cluster_datasets/PXD021518_results/result.csv'
output_file = '../data/PXD021518/maracluster/MaRaCluster_processed.clusters_p5_enriched.tsv'
# Call the function with the specified file paths
cluster_results = process_maracluster_results(input_file)
cluster_results = enrich_cluster_results(cluster_results,input_workflow_reuslt)
print(cluster_results.head())
cluster_results.to_csv(output_file, sep='\t', index=False)
