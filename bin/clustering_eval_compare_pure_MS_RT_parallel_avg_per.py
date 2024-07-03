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
from multiprocessing import Pool
from functools import partial

method_dic = {'mscluster':{'filename':'#Filename','scan':'#Scan','mass':'#ParentMass','rt_time':'#RetTime','cluster_id':'#ClusterIdx'},
              'falcon':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time','cluster_id':'cluster'},
              'maracluster':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time','cluster_id':'cluster'}}

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
                rt = float(spectrum.scan_time_in_minutes()) * 60
                mass =precursors[0]['mz']
                scan_number = int(spectrum['id'])
                ms1_data.append((mass, rt, scan_number))
    except Exception as e:
            print(e)

    return ms1_data

def are_matching_spectra(filename, scan1, scan2,matching_pairs_set):
    return (filename, scan1, scan2) in matching_pairs_set

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
           and abs(spec1[1] - spec2[1]) <= 0.01 and abs(spec1[2] - spec2[2]) <= 0.3
    ]
    G.add_edges_from(edges)

    return G
def create_matching_network(cluster, method):
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
            if abs(spec1[method_dic[method]['mass']]-spec2[method_dic[method]['mass']])<=0.01 and abs(spec1[method_dic[method]['rt_time']]-spec2[method_dic[method]['rt_time']])<=0.1:
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
def calculate_cluster_purity(cluster,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, method)

    max_component_sizes = calculate_max_component_per_file(G)
    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster[method_dic[method]['filename']])

    max_fraction = 0.0

    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        fraction = max_component_sizes[filename] /count
        max_fraction = max(fraction, max_fraction)

    return max_fraction

def calculate_cluster_purity_weighted_avg(task,method):
    _,cluster = task
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = optimized_create_matching_network(cluster, method)

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

def calculate_cluster_purity_avg(cluster,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster,method)

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
        return '1 to 1'  # Change here to ensure sorting works
    else:
        upper_limit = 2
        while size > upper_limit:
            upper_limit *= 2
        return f"{upper_limit//2+1} to {upper_limit}"


def apply_parallel_with_tqdm(groups, func, method):
    tasks = [(name, group) for name, group in groups]  # Ensure tasks are prepared as tuples
    partial_func = partial(func, method=method)

    with Pool(processes=64) as pool:
        result_list = []
        for result in tqdm(pool.imap(partial_func, tasks), total=len(tasks)):
            result_list.append(result)

    return result_list
def mscluster_purity(cluster_results):
    #handle the new version workflow filename issue
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzml')

    #return cluster_results.groupby('#ClusterIdx').progress_apply(lambda x: calculate_cluster_purity_weighted_avg(x, matching_pairs_set,'mscluster')), cluster_results.groupby('#ClusterIdx').size()
    return apply_parallel_with_tqdm(cluster_results.groupby('#ClusterIdx'),calculate_cluster_purity_weighted_avg,'mscluster'), cluster_results.groupby('#ClusterIdx').size()

def falcon_purity(cluster_results):
    #return cluster_results.groupby('cluster').progress_apply(lambda x: calculate_cluster_purity_weighted_avg(x,'falcon')), cluster_results.groupby('cluster').size()
    return apply_parallel_with_tqdm(cluster_results.groupby('cluster'), calculate_cluster_purity_weighted_avg, 'falcon'), cluster_results.groupby('cluster').size()
def maracluster_purity(cluster_results):
    #return cluster_results.groupby('cluster').progress_apply(lambda x: calculate_cluster_purity_weighted_avg(x, 'maracluster')), cluster_results.groupby('cluster').size()
    return apply_parallel_with_tqdm(cluster_results.groupby('cluster'), calculate_cluster_purity_weighted_avg,
                                    'maracluster'), cluster_results.groupby('cluster').size()
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

def mscluster_merge(cluster_results):
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzml')
    cluster_results['#RetTime'] = cluster_results['#RetTime']

    # Create a merged dataset with left join to include all cluster results and match with database results

    return cluster_results
def falcon_merge(cluster_results):

    return cluster_results
def maracluster_merge(cluster_results,reference_data):
    # reference_data.drop('#ClusterIdx', axis=1, inplace=True)
    # reference_data['#Filename'] = reference_data['#Filename'].str.replace('input_spectra', 'mzml')
    # reference_data['#RetTime'] = reference_data['#RetTime']
    # cluster_results['filename'] = cluster_results['filename'].str.replace('data','mzml')
    # merged_data = pd.merge(cluster_results, reference_data, left_on=['filename', 'scan'], right_on=['#Filename', '#Scan'],how='inner')
    return cluster_results
def custom_sort(bin_range):
    if bin_range == '1 to 1':
        return 0, 1  # This ensures "1 to 1" has the smallest possible sort key
    else:
        # Extract the start and end numbers from the bin label and use them as sort keys
        start, end = map(int, bin_range.split(' to '))
        return start, end

def range_avg(cluster_purity, cluster_size):
    purity_by_cluster_size = pd.DataFrame({'ClusterSize': cluster_size, 'ClusterPurity': cluster_purity})
    purity_by_cluster_size['Bin'] = purity_by_cluster_size['ClusterSize'].apply(assign_size_range_bin)

    # Group by 'Bin' and calculate the mean purity for each bin
    average_purity_by_bin = purity_by_cluster_size.groupby('Bin')['ClusterPurity'].mean().reset_index()

    count_clusters_by_bin = purity_by_cluster_size.groupby('Bin')['ClusterSize'].count().reset_index()
    count_clusters_by_bin.rename(columns={'ClusterSize': 'ClusterCount'}, inplace=True)

    # Merge the two dataframes on 'Bin'
    merged_df = pd.merge(average_purity_by_bin, count_clusters_by_bin, on='Bin')

    # Sort by 'Bin'
    merged_df = merged_df.sort_values('Bin', key=lambda x: x.map(custom_sort))

    return merged_df

def calculate_weighted_average_purity(purity, cluster_size):
    weighted_purity = (purity * cluster_size).sum() / cluster_size.sum()
    return weighted_purity

def calculate_n50(cluster_size, total_spectra):
    sorted_sizes = np.sort(cluster_size)[::-1]  # Sort cluster sizes in descending order
    #total_spectra = sorted_sizes.sum()
    cumulative_sum = 0
    n50 = None
    for size in sorted_sizes:
        cumulative_sum += size
        if cumulative_sum >= total_spectra *0.1:
            n50 = size
            break
    return n50


def calculate_clustered_percentage(clustering_results: pd.DataFrame, total_spectra: int, method: str) -> float:
    """
    Calculate the percentage of clustered spectra.

    Parameters:
    clustering_results (pd.DataFrame): DataFrame containing clustering results with a column named 'Cluster_Id'.
    total_spectra (int): Total number of spectra in the dataset.

    Returns:
    float: Percentage of clustered spectra.
    """
    # Count the number of spectra in each cluster
    cluster_counts = clustering_results[method_dic[method]['cluster_id']].value_counts()

    # Filter clusters that have size >= 2
    clustered_spectra = cluster_counts[cluster_counts >= 2].sum()

    # Calculate the percentage
    percentage_clustered = (clustered_spectra / total_spectra) * 100

    return percentage_clustered


if __name__ == "__main__":
    tqdm.pandas()

    # Old data path
    # mscluster_results = pd.read_csv('../data/results/nf_output/clustering/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    # falcon_results = pd.read_csv('../data/cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    #print("original falcon_results:",len(falcon_results))
    # maracluster_results= pd.read_csv('./processed_clusters.tsv', sep='\t')

    # MSV000081981 default settings path
    # mscluster_results = pd.read_csv('../data/MSV000081981/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    # falcon_results = pd.read_csv('../data/MSV000081981/Falcon_cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    # maracluster_results= pd.read_csv('/home/user/LabData/XianghuData/metabolomics_clustering_data/MSV000081981_results/maracluster/MaRaCluster_processed.clusters_p5_enriched.tsv', sep='\t')

    mscluster_results = pd.read_csv('../data/MSV000081981/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    falcon_results = pd.read_csv('/home/user/LabData/XianghuData/Falcon_Cluster_Benchmark/MSV000081981_0.2/output_summary/cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    maracluster_results= pd.read_csv('/home/user/LabData/XianghuData/metabolomics_clustering_data/MSV000081981_results/maracluster/MaRaCluster_processed.clusters_p15_enriched.tsv', sep='\t')



    mscluster_percentage = calculate_clustered_percentage(mscluster_results, 7838866,'mscluster')
    falcon_percentage = calculate_clustered_percentage(falcon_results, 7838866,'falcon')
    # online_falcon_n50 = calculate_n50(online_falcon_size,109333)
    maracluster_percentage = calculate_clustered_percentage(maracluster_results, 7838866,'maracluster')

    print("Clustered percentage for MSCluster:", mscluster_percentage)
    print("Clustered percentage for Falcon:", falcon_percentage)
    # print("N50 value for Online Falcon:", online_falcon_n50)
    print("Clustered percentage for MaRaCluster:", maracluster_percentage)

    # database_results = pd.read_csv('./PXD021518_filtered.tsv', sep='\t')  # Adjust file path and format accordingly
    data_folder_path = os.path.dirname(os.getcwd()) + '/data'
    matching_pairs_set_filename = 'matching_pairs_set.pkl'
    mscluster_results_copy = mscluster_results.copy()
    mscluster_results = mscluster_merge(mscluster_results)
    falcon_results = falcon_merge(falcon_results)
    print("mscluster cluster number:", len(mscluster_results.groupby('#ClusterIdx').size()))
    cluster_sizes = mscluster_results.groupby('#ClusterIdx').size()

    # Filter clusters larger than 4096
    clusters_larger_than_4096 = cluster_sizes[cluster_sizes > 4097]

    print(f"Number of mscluster clusters larger than 4097: {len(clusters_larger_than_4096)}")
    print("falcon cluster number:", len(falcon_results.groupby('cluster').size()))
    print("maracluster cluster number:", len(maracluster_results.groupby('cluster').size()))
    cluster_sizes = maracluster_results.groupby('cluster').size()

    # Filter clusters larger than 4096
    clusters_larger_than_4096 = cluster_sizes[cluster_sizes > 4097]

    print(f"Number of mscluster clusters larger than 4097: {len(clusters_larger_than_4096)}")
    maracluster_results = maracluster_merge(maracluster_results,mscluster_results_copy)
    print("start calculating the mscluster results")
    mscluster_purity, mscluster_size = mscluster_purity(mscluster_results)
    print("start calculating the Falcon results")
    falcon_purity, falcon_size = falcon_purity(falcon_results)
    print("start calculating the Maracluster results")
    maracluster_purity, maracluster_size = maracluster_purity(maracluster_results)
    print('falcon size:',len(falcon_results))
    print('mscluster size:',len(mscluster_results))
    print('maracluster size:',len(maracluster_results))

    mscluster_weighted_avg_purity = sum(mscluster_purity)/len(mscluster_purity)
    falcon_weighted_avg_purity = sum(falcon_purity)/len(falcon_purity)
    # online_falcon_weighted_avg_purity = calculate_weighted_average_purity(online_falcon_purity, online_falcon_size)
    maracluster_weighted_avg_purity = sum(maracluster_purity)/len(maracluster_purity)

    # mscluster_weighted_avg_purity = calculate_weighted_average_purity(mscluster_purity, mscluster_size)
    # falcon_weighted_avg_purity = calculate_weighted_average_purity(falcon_purity, falcon_size)
    # # online_falcon_weighted_avg_purity = calculate_weighted_average_purity(online_falcon_purity, online_falcon_size)
    # maracluster_weighted_avg_purity = calculate_weighted_average_purity(maracluster_purity, maracluster_size)

    print("Weighted Average Purity for MSCluster:", mscluster_weighted_avg_purity)
    print("Weighted Average Purity for Falcon:", falcon_weighted_avg_purity)
    # print("Weighted Average Purity for Online Falcon:", online_falcon_weighted_avg_purity)
    print("Weighted Average Purity for MaRaCluster:", maracluster_weighted_avg_purity)

    # methods = ['MSCluster', 'Falcon', 'Online Falcon', 'MaRaCluster']
    # n50_values = [mscluster_n50, falcon_n50, online_falcon_n50, maracluster_n50]
    # weighted_avg_purities = [mscluster_weighted_avg_purity, falcon_weighted_avg_purity, online_falcon_weighted_avg_purity, maracluster_weighted_avg_purity]



