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
           and abs(spec1[1] - spec2[1]) <= 0.01 and abs(spec1[2] - spec2[2]) <= 0.49
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
            if abs(spec1[method_dic[method]['mass']]-spec2[method_dic[method]['mass']])<=0.01 and abs(spec1[method_dic[method]['rt_time']]-spec2[method_dic[method]['rt_time']])<=0.5:
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

if __name__ == "__main__":
    folder_path = '/home/user/LabData/XianghuData/MS_Cluster_datasets/Combine_test/mzML'
    mscluster_results = pd.read_csv('../data/results/nf_output/clustering/clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    falcon_results = pd.read_csv('../data/PXD023047/falcon/falcon_cluster_info_0.3.tsv',sep='\t')  # Adjust file path and format accordingly
    maracluster_results= pd.read_csv('../data/PXD023047/maracluster/MaRaCluster_processed.clusters_p5_enriched.tsv', sep='\t')
    # falcon_results = pd.read_csv('../data/cluster_info.tsv', sep='\t')  # Adjust file path and format accordingly
    # maracluster_results= pd.read_csv('../data/PXD023047/maracluster/MaRaCluster_processed.clusters_p10_enriched.tsv', sep='\t')
    mscluster_n50 = calculate_n50(mscluster_results.groupby('#ClusterIdx').size(), 109333)
    falcon_n50 = calculate_n50(falcon_results.groupby('cluster').size(), 109333)
    # online_falcon_n50 = calculate_n50(online_falcon_size,109333)
    maracluster_n50 = calculate_n50(maracluster_results.groupby('cluster').size(), 109333)

    # mscluster_results = pd.read_csv('../data/Combine_results/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    # falcon_results = pd.read_csv('../data/Combine_results/Falcon_cluster_info.tsv',sep='\t')  # Adjust file path and format accordingly
    # maracluster_results= pd.read_csv('../data/Combine_results/MaRaCluster_processed.clusters_p10_enriched.tsv', sep='\t')
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
    mscluster_weighted_avg_purity = calculate_weighted_average_purity(mscluster_purity, mscluster_size)
    falcon_weighted_avg_purity = calculate_weighted_average_purity(falcon_purity, falcon_size)
    # online_falcon_weighted_avg_purity = calculate_weighted_average_purity(online_falcon_purity, online_falcon_size)
    maracluster_weighted_avg_purity = calculate_weighted_average_purity(maracluster_purity, maracluster_size)

    print("Weighted Average Purity for MSCluster:", mscluster_weighted_avg_purity)
    print("Weighted Average Purity for Falcon:", falcon_weighted_avg_purity)
    # print("Weighted Average Purity for Online Falcon:", online_falcon_weighted_avg_purity)
    print("Weighted Average Purity for MaRaCluster:", maracluster_weighted_avg_purity)


    print("N50 value for MSCluster:", mscluster_n50)
    print("N50 value for Falcon:", falcon_n50)
    # print("N50 value for Online Falcon:", online_falcon_n50)
    print("N50 value for MaRaCluster:", maracluster_n50)

    # methods = ['MSCluster', 'Falcon', 'Online Falcon', 'MaRaCluster']
    # n50_values = [mscluster_n50, falcon_n50, online_falcon_n50, maracluster_n50]
    # weighted_avg_purities = [mscluster_weighted_avg_purity, falcon_weighted_avg_purity, online_falcon_weighted_avg_purity, maracluster_weighted_avg_purity]
    methods = ['MSCluster', 'Falcon', 'MaRaCluster']
    n50_values = [mscluster_n50, falcon_n50, maracluster_n50]
    weighted_avg_purities = [mscluster_weighted_avg_purity, falcon_weighted_avg_purity, maracluster_weighted_avg_purity]

    # falcon_x= [4,17,18,18,18]
    falcon_x = [8, 24, 25, 26, 26]
    falcon_y = [0.9896683397928712,0.9748935149079129,0.973328106397189,0.972145244693312,0.9704464146834656]

    # maracluster_x = [17,17,16,14,11]
    maracluster_x = [25,23, 21, 18, 16, 14]
    maracluster_y = [0.972550831792976,0.9771322946923686,0.9806311064166887,0.9839054660681278,0.9875231053604436,0.9911407446527595]

    falcon_x_DB = [8, 24, 25, 26, 26]
    falcon_y_DB = [0.995044158, 0.98821705, 0.986540639, 0.984972717, 0.983399087]

    # maracluster_x = [17,17,16,14,11]
    maracluster_x_DB = [25,23, 21, 18, 16, 14]
    maracluster_y_DB = [0.9845260100343279, 0.98844732, 0.990071297, 0.991444415, 0.992632691, 0.993992606]

    mscluster_x_DB = [40]
    mscluster_y_DB = [0.9362784845519363]


    # Plot
    # plt.figure(figsize=(10, 6))
    # valid_points = [(n50, purity) for purity, n50 in zip(weighted_avg_purities, n50_values) if
    #                 purity is not None and n50 is not None]
    # valid_methods = [method for method, purity, n50 in zip(methods, weighted_avg_purities, n50_values) if
    #                  purity is not None and n50 is not None]
    # valid_n50, valid_purities = zip(*valid_points)
    # plt.scatter(valid_n50, valid_purities, color='blue', label='Methods_Default_MSRT')
    # plt.scatter(falcon_x, falcon_y, color='red', label='Falcon_MSRT')
    # plt.scatter(maracluster_x, maracluster_y, color='green', label='MaraCluster_MSRT')
    # plt.scatter(falcon_x_DB, falcon_y_DB, label='Falcon_DB')
    # plt.scatter(maracluster_x_DB, maracluster_y_DB, label='MaraCluster_DB')
    # plt.scatter(mscluster_x_DB, mscluster_y_DB, label='MSCluster_DB')
    # plt.xlabel('N50 Value',fontsize=12)
    # plt.ylabel('Weighted Average Purity',fontsize=12)
    # #plt.title('N50 vs. Weighted Average Purity for Each Method')
    # plt.ylim([0.85,1.01])
    # for i, method in enumerate(valid_methods):
    #     plt.annotate(method, (valid_n50[i], valid_purities[i]), textcoords="offset points", xytext=(0, 10), ha='center',fontsize=12)
    # plt.legend(loc='lower left',fontsize=12)
    # plt.grid(True)
    # plt.tight_layout()
    # plt.show()
    plt.figure(figsize=(10, 6))

    # Create valid points
    valid_points = [(purity, n50) for purity, n50 in zip(weighted_avg_purities, n50_values) if
                    purity is not None and n50 is not None]
    valid_methods = [method for method, purity, n50 in zip(methods, weighted_avg_purities, n50_values) if
                     purity is not None and n50 is not None]

    # Unpack valid points into two separate lists
    valid_purities, valid_n50 = zip(*valid_points)

    # Plot scatter points
    #plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
    plt.plot(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
    plt.plot(maracluster_y, maracluster_x,'o-', color='green', label='MaraCluster_MSRT')
    plt.plot(falcon_y_DB, falcon_x_DB,'o-', label='Falcon_DB')
    plt.plot(maracluster_y_DB, maracluster_x_DB,'o-', label='MaraCluster_DB')
    plt.scatter(mscluster_y_DB, mscluster_x_DB, label='MSCluster_DB')
    plt.scatter(0.9283237802320704,40, label='MSCluster_MSRT')


    # Adjust axis labels
    plt.xlabel('Weighted Average Purity', fontsize=12)
    plt.ylabel('N50 Value', fontsize=12)

    # Reverse the x-axis
    plt.gca().invert_xaxis()


    # Annotate points with method names
    # for i, method in enumerate(valid_methods):
    #     plt.annotate(method, (valid_purities[i], valid_n50[i]), textcoords="offset points", xytext=(0, 10), ha='center',
    #                  fontsize=12)

    # Add legend, grid, and layout adjustments
    plt.legend(loc='lower right', fontsize=12)
    plt.title("Clustering Benchmark Results or Proteomics Dataset")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
