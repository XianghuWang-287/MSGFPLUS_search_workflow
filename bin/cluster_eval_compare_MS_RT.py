import os
from itertools import combinations
from pyteomics import mzml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymzml
from  xml.etree.ElementTree import ParseError
import networkx as nx
from collections import Counter
from matplotlib.colors import LogNorm
import pickle
import copy
from tqdm import tqdm
from collections import defaultdict

method_dic = {'mscluster':{'filename':'#Filename','scan':'#Scan','mass':'#ParentMass','rt_time':'#RetTime'},
              'falcon':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time'},
              'maracluster':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time'}}

def get_ms1_data(mzml_file):
    ms1_data = []

    run = pymzml.run.Reader(mzml_file, build_index_from_scratch=False)



    try:
        for idx, spectrum in enumerate(run):
            if spectrum.ms_level == 2:
                precursors = spectrum.selected_precursors
                if len(precursors) != 1:
                    exception = NotImplementedError(
                        "Expected one precursors but got {} while checking index {}.".format(len(precursors),
                                                                                             spectrum.index))
                    print(exception)
                    print(mzml_file)
                    continue
                rt = float(spectrum.scan_time_in_minutes())
                mass =precursors[0]['mz']
                scan_number = int(spectrum['id'])
                ms1_data.append((mass, rt, scan_number))
    except Exception as e:
            print(e)

    return ms1_data

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

def are_matching_spectra(filename, scan1, scan2,matching_pairs_set):
    return (filename, scan1, scan2) in matching_pairs_set

from itertools import combinations
from collections import defaultdict

def optimized_create_matching_network(cluster, method):
    G = nx.Graph()

    # Precompute the node names and add them to the graph
    node_attrs = {
        f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}": {"filename": row[method_dic[method]['filename']]}
        for _, row in cluster.iterrows()
    }
    G.add_nodes_from(node_attrs.items())

    # Precompute mass and rt_time for efficient access
    specs = cluster.apply(lambda row: (
        f"{row[method_dic[method]['filename']]}_{row[method_dic[method]['scan']]}",
        row[method_dic[method]['mass']],
        row[method_dic[method]['rt_time']]
    ), axis=1).tolist()

    # Create edges based on conditions
    edges = [
        (spec1[0], spec2[0]) for spec1, spec2 in combinations(specs, 2)
        if spec1[0].split('_')[0] == spec2[0].split('_')[0]  # Ensure filenames are the same
           and abs(spec1[1] - spec2[1]) <= 0.01 and abs(spec1[2] - spec2[2]) <= 12
    ]
    G.add_edges_from(edges)

    return G

def create_matching_network(cluster, matching_pairs_set,method):
    G = nx.Graph()

    row_indices = range(len(cluster))
    row_combinations = combinations(row_indices, 2)
    for _, spec in cluster.iterrows():
        node_a = f"{spec[method_dic[method]['filename']]}_{spec[method_dic[method]['scan']]}"
        G.add_node(node_a, filename=spec[method_dic[method]['filename']])
    for combo in row_combinations:
        index1, index2 = combo
        spec1 = cluster.iloc[index1]
        spec2 = cluster.iloc[index2]
        if spec1[method_dic[method]['filename']] == spec2[method_dic[method]['filename']]:
            if abs(spec1[method_dic[method]['mass']]-spec2[method_dic[method]['mass']])<=0.01 and abs(spec1[method_dic[method]['rt_time']]-spec2[method_dic[method]['rt_time']])<=10:
                node_a = f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
                node_b = f"{spec2[method_dic[method]['filename']]}_{spec2[method_dic[method]['scan']]}"
                G.add_edge(node_a,node_b)

    return G

def create_matching_network_across_files(cluster,method,tolerence):
    G = nx.Graph()

    for _, spec1 in cluster.iterrows():
        node_a =f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
        G.add_node(node_a, filename=spec1[method_dic[method]['filename']])
        for _, spec2 in cluster.iterrows():
            if (spec1[method_dic[method]['filename']] == spec2[method_dic[method]['filename']] and spec1[method_dic[method]['scan']] == spec2[method_dic[method]['scan']]):
                continue
            if abs(spec1['']):
                node_a =f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
                node_b =f"{spec2[method_dic[method]['filename']]}_{spec2[method_dic[method]['scan']]}"
                G.add_node(node_a,filename=spec1[method_dic[method]['filename']])
                G.add_node(node_b,filename=spec2[method_dic[method]['filename']])
                G.add_edge(node_a,node_b)
    return G


def calculate_max_component_per_file(G):

    # Find all connected components in the graph
    components = nx.connected_components(G)

    # Initialize a dictionary to hold the maximum component size for each file
    max_component_sizes = defaultdict(int)

    # Iterate through each component
    for component in components:
        # Create a temporary dictionary to count the number of nodes per file in this component
        file_counts = defaultdict(int)

        # Count nodes per file in the current component
        for node in component:
            filename = G.nodes[node]['filename']
            file_counts[filename] += 1

        # Update the max component size for each file encountered in this component
        for filename, count in file_counts.items():
            if count > max_component_sizes[filename]:
                max_component_sizes[filename] = count
    # Ensure that files represented by single nodes are accounted for
    for node in G.nodes:
        filename = G.nodes[node]['filename']
        if filename not in max_component_sizes:
            max_component_sizes[filename] = 1
        else:
            # Ensure there's at least a count of 1 for each file,
            # useful if a file's node(s) were not part of any component counted above
            max_component_sizes[filename] = max(max_component_sizes[filename], 1)

    return dict(max_component_sizes)

# Function to calculate purity for a cluster
def calculate_cluster_purity(cluster, matching_pairs_set,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = optimized_create_matching_network(cluster,method)

    max_component_sizes = calculate_max_component_per_file(G)
    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    max_fraction = 0.0

    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        fraction = max_component_sizes[filename] /count
        max_fraction = max(fraction, max_fraction)

    return max_fraction

def calculate_cluster_purity_weighted_avg(cluster, matching_pairs_set,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = optimized_create_matching_network(cluster,method)

    max_component_sizes = calculate_max_component_per_file(G)

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max_component_sizes[filename]
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    weighted_sum = sum(value * frequency for value, frequency in zip(values, frequencies))

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    weighted_average = weighted_sum / total_frequency
    return weighted_average

def calculate_cluster_purity_avg(cluster, matching_pairs_set,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, matching_pairs_set,method)

    # Get connected components
    max_component_sizes = calculate_max_component_per_file(G)
    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max_component_sizes[filename]
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    total_purity = sum(values)

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    average = total_purity / total_frequency
    return average


def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=6000):
    mass_diff = abs(scan1[0] - scan2[0])
    rt_diff = abs(scan1[1] - scan2[1])
    return mass_diff <= mass_tolerance and rt_diff <= rt_tolerance

def assign_size_range_bin(size):
    if size == 1:
        return '1'
    else:
        upper_limit = 2
        while size > upper_limit:
            upper_limit *= 2
        return f"{upper_limit//2+1} to {upper_limit}"

def mscluster_purity(cluster_results,matching_pairs_set):
    #handle the new version workflow filename issue
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzML')

    return cluster_results.groupby('#ClusterIdx').apply(lambda x: calculate_cluster_purity_weighted_avg(x, matching_pairs_set,'mscluster')), cluster_results.groupby('#ClusterIdx').size()

def falcon_purity(cluster_results,matching_pairs_set):
    return cluster_results.groupby('cluster').apply(lambda x: calculate_cluster_purity_weighted_avg(x, matching_pairs_set,'falcon')), cluster_results.groupby('cluster').size()
def maracluster_purity(cluster_results,matching_pairs_set):
    return cluster_results.groupby('cluster').apply(
        lambda x: calculate_cluster_purity_weighted_avg(x, matching_pairs_set, 'maracluster')), cluster_results.groupby(
        'cluster').size()
def generate_matching_pairs_set_file(folder_path,data_folder_path):
    matching_pairs_all_files = []
    for filename in tqdm(os.listdir(folder_path)):
        if filename.endswith('.mzML'):
            mzml_file = os.path.join(folder_path, filename)
            ms1_data = get_ms1_data(mzml_file)

            matching_pairs = []

            for scan1, scan2 in combinations(ms1_data, 2):
                if compare_scans(scan1, scan2):
                    matching_pairs.append((scan1[2], scan2[2]))

            matching_pairs_all_files.extend([(filename, pair[0], pair[1]) for pair in matching_pairs])

    matching_pairs_all_files = [(f'mzML/{filename}', scan_start, scan_end) for filename, scan_start, scan_end in
                                matching_pairs_all_files]
    matching_pairs_set = {(item[0], item[1], item[2]) for item in matching_pairs_all_files}
    with open(data_folder_path+'/matching_pairs_set.pkl', 'wb') as file:
        pickle.dump(matching_pairs_set, file)
    return matching_pairs_set

def mscluster_merge(cluster_results,database_results):
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzML')
    cluster_results['#RetTime'] = cluster_results['#RetTime']
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
    # Create a merged dataset with left join to include all cluster results and match with database results
    merged_data = cluster_results.merge(database_results, left_on=['#Filename', '#Scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    return merged_data
def falcon_merge(cluster_results,database_results):
    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    merged_data.to_csv("falcon_merge.tsv", sep='\t', index=False)
    return merged_data
def maracluster_merge(cluster_results,reference_data,database_results):
    reference_data.drop('#ClusterIdx', axis=1, inplace=True)
    reference_data['#Filename'] = reference_data['#Filename'].str.replace('input_spectra', 'mzML')
    reference_data['#RetTime'] = reference_data['#RetTime']
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
    merged_data = pd.merge(cluster_results, reference_data, left_on=['filename', 'scan'], right_on=['#Filename', '#Scan'],how='inner')
    merged_data = merged_data.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    return merged_data

def maracluster_merge_new(cluster_results,database_results):
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']

    merged_data = cluster_results.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    return merged_data


if __name__ == "__main__":
    folder_path = '/home/user/LabData/XianghuData/MS_Cluster_datasets/Combine_test/mzML'
    # mscluster_results = pd.read_csv('../data/results/nf_output/clustering/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    # falcon_results = pd.read_csv('../data/cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    # maracluster_results= pd.read_csv('./processed_clusters.tsv', sep='\t')

    mscluster_results = pd.read_csv('../data/Combine_results/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    falcon_results = pd.read_csv('../data/Combine_results/Falcon_cluster_info.tsv',sep='\t')  # Adjust file path and format accordingly
    maracluster_results= pd.read_csv('../data/Combine_results/MaRaCluster_processed.clusters_p10_enriched.tsv', sep='\t')
    print(maracluster_results)
    print("original mscluster spectra number:", len(mscluster_results))
    print("original falcon spectra number:", len(falcon_results))
    print("original maracluster spectra number:", len(maracluster_results))
    print("mscluster cluster number before merge:", len(mscluster_results.groupby('#ClusterIdx').size()))
    print("falcon cluster number before merge:", len(falcon_results.groupby('cluster').size()))
    print("maracluster cluster number before merge:", len(maracluster_results.groupby('cluster').size()))
    database_results = pd.read_csv('./Combine_test_filtered.tsv', sep='\t')  # Adjust file path and format accordingly
    #database_results = database_results[database_results['DB:EValue'] < 0.002]
    database_results['Peptide'] = database_results['Peptide'].str.replace('I', 'L')
    print("number of different piptides: ",len(database_results.groupby('Peptide').size()))


    data_folder_path = os.path.dirname(os.getcwd()) + '/data'
    matching_pairs_set_filename = 'matching_pairs_set.pkl'
    if matching_pairs_set_filename not in os.listdir(data_folder_path):
        print('matching_pairs_set_filename not exist.')
        print('Creating matching_pairs_set_filename...')
        matching_pairs_set = generate_matching_pairs_set_file(folder_path,data_folder_path)
    else:
        with open(data_folder_path+'/'+matching_pairs_set_filename, 'rb') as file:
            matching_pairs_set = pickle.load(file)
    mscluster_results_copy = mscluster_results.copy()
    mscluster_results = mscluster_merge(mscluster_results,copy.copy(database_results))
    falcon_results = falcon_merge(falcon_results,copy.copy(database_results))
    maracluster_results = maracluster_merge_new(maracluster_results,copy.copy(database_results))
    print("merged mscluster spectra number:", len(mscluster_results))
    print("merged falcon spectra number:", len(falcon_results))
    print("merged maracluster spectra number:", len(maracluster_results))
    print("mscluster cluster number after merge:", len(mscluster_results.groupby('#ClusterIdx').size()))
    print("falcon cluster number after merge:", len(falcon_results.groupby('cluster').size()))
    print("maracluster cluster number after merge:",len(maracluster_results.groupby('cluster').size()))
    mscluster_purity, mscluster_size = mscluster_purity(mscluster_results, matching_pairs_set)
    falcon_purity, falcon_size = falcon_purity(falcon_results,matching_pairs_set)
    maracluster_purity, maracluster_size = maracluster_purity(maracluster_results,matching_pairs_set)
    print('falcon size:',len(falcon_results))
    print('mscluster size:',len(mscluster_results))
    print('maracluster size:',len(maracluster_results))
    purity_by_cluster_size_mscluster = pd.DataFrame({'ClusterSize': mscluster_size, 'ClusterPurity': mscluster_purity})
    purity_by_cluster_size_mscluster['Bin'] = purity_by_cluster_size_mscluster['ClusterSize'].apply(
        assign_size_range_bin)
    # Group by cluster size and calculate average purity for each group
    average_purity_by_bin_mscluster = purity_by_cluster_size_mscluster.groupby('Bin')[
        'ClusterPurity'].mean().reset_index()

    # Since the bins are non-numeric and custom, we'll sort them by the lower bound of each range for plotting
    average_purity_by_bin_mscluster['SortKey'] = average_purity_by_bin_mscluster['Bin'].apply(
        lambda x: int(x.split(' ')[0]))
    average_purity_by_bin_mscluster = average_purity_by_bin_mscluster.sort_values('SortKey')

    purity_by_cluster_size_falcon = pd.DataFrame({'ClusterSize': falcon_size, 'ClusterPurity': falcon_purity})
    purity_by_cluster_size_falcon['Bin'] = purity_by_cluster_size_falcon['ClusterSize'].apply(assign_size_range_bin)
    average_purity_by_bin_falcon = purity_by_cluster_size_falcon.groupby('Bin')['ClusterPurity'].mean().reset_index()
    average_purity_by_bin_falcon['SortKey'] = average_purity_by_bin_falcon['Bin'].apply(lambda x: int(x.split(' ')[0]))
    average_purity_by_bin_falcon = average_purity_by_bin_falcon.sort_values('SortKey')

    purity_by_cluster_size_maracluster = pd.DataFrame({'ClusterSize': maracluster_size, 'ClusterPurity': maracluster_purity})
    purity_by_cluster_size_maracluster['Bin'] = purity_by_cluster_size_maracluster['ClusterSize'].apply(assign_size_range_bin)
    average_purity_by_bin_maracluster = purity_by_cluster_size_maracluster.groupby('Bin')['ClusterPurity'].mean().reset_index()
    average_purity_by_bin_maracluster['SortKey'] = average_purity_by_bin_maracluster['Bin'].apply(lambda x: int(x.split(' ')[0]))
    average_purity_by_bin_maracluster =average_purity_by_bin_maracluster.sort_values('SortKey')
    print(average_purity_by_bin_mscluster)
    print(average_purity_by_bin_falcon)
    print(average_purity_by_bin_maracluster)
    # Plotting
    plt.figure(figsize=(12, 6))
    plt.plot(average_purity_by_bin_mscluster['Bin'], average_purity_by_bin_mscluster['ClusterPurity'], marker='o',
             linestyle='-', label='mscluster')
    plt.plot(average_purity_by_bin_falcon['Bin'], average_purity_by_bin_falcon['ClusterPurity'], marker='o',
             linestyle='-', label='falcon')
    plt.plot(average_purity_by_bin_maracluster['Bin'],average_purity_by_bin_maracluster['ClusterPurity'], marker='o',linestyle='-', label='maracluster')
    plt.xlabel('Cluster Size Range')
    plt.ylabel('Average Purity')
    plt.title('Average Cluster Purity by Size Range MS-RT')
    plt.xticks(rotation=45)
    plt.ylim([0,1.1])
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.legend()
    plt.tight_layout()  # Adjust layout to make room for the rotated x-axis labels
    plt.show()
    plt.savefig('combine_test.png')

