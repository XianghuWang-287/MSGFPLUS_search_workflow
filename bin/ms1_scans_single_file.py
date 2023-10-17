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

# def get_ms1_data(mzml_file):
#     ms1_data = []
#     with mzml.read(mzml_file) as reader:
#         for spectrum in reader:
#             if spectrum['ms level'] == 2:
#                 mass = float(spectrum['precursorList']['precursor']['selectedIonList']['cvParam'].get('value', 0))
#                 rt = float(spectrum['scanList']['scan'][0]['scan start time'])
#                 scan_number = int(spectrum['id'].split('scan=')[1])
#                 ms1_data.append((mass, rt,scan_number))
#     return ms1_data

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


def create_matching_network(cluster, matching_pairs_set):
    G = nx.Graph()

    for _, spec1 in cluster.iterrows():
        for _, spec2 in cluster.iterrows():
            if spec1['#Filename'] != spec2['#Filename']:
                continue

            if are_matching_spectra(spec1['#Filename'], spec1['#Scan'], spec2['#Scan'], matching_pairs_set):
                node_a =f"{spec1['#Filename']}_{spec1['#Filename']}"
                node_b =f"{spec2['#Filename']}_{spec2['#Filename']}"
                G.add_node(node_a,filename=spec1['#Filename'])
                G.add_node(node_b,filename=spec2['#Filename'])
                G.add_edge(node_a,node_b)

    return G

# Function to calculate purity for a cluster
def calculate_cluster_purity(cluster, matching_pairs_set):
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, matching_pairs_set)

    # Get connected components
    components = list(nx.connected_components(G))

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster['#Filename'])

    max_fraction = 0.0

    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max((comp.number_of_nodes() for comp in S if any(comp.nodes[node]['filename'] == filename for node in comp.nodes())), default=1)
        fraction = largest_component_size / count
        if fraction>1:
            print("total_size",G.number_of_nodes())
            print("largeest component size",largest_component_size)
        max_fraction = max(max_fraction, fraction)

    return max_fraction



def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=30):
    mass_diff = abs(scan1[0] - scan2[0])
    rt_diff = abs(scan1[1] - scan2[1])
    return mass_diff <= mass_tolerance and rt_diff <= rt_tolerance


if __name__ == "__main__":
    folder_path = '/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_convert/mzML'

    cluster_results = pd.read_csv('/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_results/msculster_results/clustering/clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly
    matching_pairs_all_files = []

    for filename in os.listdir(folder_path):
        if filename.endswith('.mzML'):
            mzml_file = os.path.join(folder_path, filename)
            ms1_data = get_ms1_data(mzml_file)

            matching_pairs = []

            for scan1, scan2 in combinations(ms1_data, 2):
                if compare_scans(scan1, scan2):
                    matching_pairs.append((scan1[2], scan2[2]))

            matching_pairs_all_files.extend([(filename, pair[0], pair[1]) for pair in matching_pairs])

    matching_pairs_all_files = [(f'mzML/{filename}', scan_start, scan_end) for filename, scan_start, scan_end in matching_pairs_all_files]

    matching_pairs_set = {(item[0], item[1], item[2]) for item in matching_pairs_all_files}

    cluster_purity = cluster_results.groupby('#ClusterIdx').apply(lambda x: calculate_cluster_purity(x, matching_pairs_set))

    cluster_size = cluster_results.groupby('#ClusterIdx').size()

    average_purity_by_cluster_size = cluster_purity.groupby(cluster_size).mean()

    plt.figure(figsize=(10, 6))
    plt.plot(average_purity_by_cluster_size.index, average_purity_by_cluster_size.values, marker='o')
    plt.xlabel('Cluster Size (Number of Spectra)')
    plt.ylabel('Average Cluster Purity')
    plt.title('Average Cluster Purity by Cluster Size')
    plt.grid(True)
    plt.show()
