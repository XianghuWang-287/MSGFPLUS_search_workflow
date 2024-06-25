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

def extract_ms1_data(mzml_file):
    ms1_data = {}
    run = pymzml.run.Reader(mzml_file, build_index_from_scratch=False)
    for spectrum in run:
        if spectrum.ms_level == 2:  # Assuming you are interested in MS2 scans for MaRaCluster results
            scan_number = int(spectrum['id'])
            rt = float(spectrum.scan_time_in_minutes())
            precursors = spectrum.selected_precursors
            if precursors:
                mass = precursors[0]['mz']
                ms1_data[scan_number] = (mass, rt)
    return ms1_data

def enrich_cluster_results(cluster_results,mzml_folder_path):
    cluster_results['precursor_mz'] = pd.NA
    cluster_results['retention_time'] = pd.NA
    for filename in os.listdir(mzml_folder_path):
        if filename.endswith(".mzML"):
            mzml_file_path = os.path.join(mzml_folder_path, filename)
            ms1_data = extract_ms1_data(mzml_file_path)
            for index, row in cluster_results.iterrows():
                if row['filename'] == "mzML/"+filename:
                    scan_number = row['scan']
                    if scan_number in ms1_data:
                        cluster_results.at[index, 'precursor_mz'] = ms1_data[scan_number][0]
                        cluster_results.at[index, 'retention_time'] = ms1_data[scan_number][1]
    return cluster_results

# Specify the input and output file paths
input_file = '/home/user/LabData/XianghuData/MS_Cluster_datasets/Combine_test/maracluster_output/MaRaCluster.clusters_p5.tsv'
output_file = '../data/Combine_results/maracluster/MaRaCluster_processed.clusters_p5_enriched.tsv'
folder_path = '/home/user/LabData/XianghuData/MS_Cluster_datasets/Combine_test/mzML'
# Call the function with the specified file paths
cluster_results = process_maracluster_results(input_file)
print(cluster_results.head())
cluster_results.to_csv(output_file, sep='\t', index=False)
