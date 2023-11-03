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
    if len(components) ==0:
        return 0

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster['#Filename'])

    max_fraction = 0.0

    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max((comp.number_of_nodes() for comp in S if any(comp.nodes[node]['filename'] == filename for node in comp.nodes())), default=1)
        fraction = largest_component_size / count
        max_fraction = max(max_fraction, fraction)

    return max_fraction

def calculate_cluster_purity_weighted_avg(cluster, matching_pairs_set):
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, matching_pairs_set)

    # Get connected components
    components = list(nx.connected_components(G))

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster['#Filename'])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max((comp.number_of_nodes() for comp in S if any(comp.nodes[node]['filename'] == filename for node in comp.nodes())), default=1)
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    weighted_sum = sum(value * frequency for value, frequency in zip(values, frequencies))

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    weighted_average = weighted_sum / total_frequency
    return weighted_average

def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=30):
    mass_diff = abs(scan1[0] - scan2[0])
    rt_diff = abs(scan1[1] - scan2[1])
    return mass_diff <= mass_tolerance and rt_diff <= rt_tolerance

if __name__ == "__main__":
    #constructing the ms1 results
    database_results = pd.read_csv('./filtered.tsv', sep='\t')  # Adjust file path and format accordingly
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

    cluster_indices = cluster_results['#ClusterIdx'].unique()

    purity_by_cluster_index = pd.DataFrame({'Cluster_indices': cluster_indices, 'ClusterPurity': cluster_purity})

    #constructing the db results
    database_results = pd.read_csv('./filtered.tsv', sep='\t')  # Adjust file path and format accordingly
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
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

    cluster_purity_db = merged_data.groupby('#ClusterIdx')['Peptide'].apply(
        lambda x: x.value_counts().max() / len(x))

    # Calculate the number of spectra in each cluster
    cluster_indices_db = merged_data['#ClusterIdx'].unique()

    purity_by_cluster_index_db = pd.DataFrame({'Cluster_indices': cluster_indices_db, 'ClusterPurity': cluster_purity_db})

    merged_purity_df = purity_by_cluster_index_db.merge(purity_by_cluster_index, on='Cluster_indices', suffixes=('_db', '_ms1'), how='inner')

    # Extract the common cluster indices and ClusterPurity values
    common_cluster_indices = merged_purity_df['Cluster_indices']
    cluster_purity_db = merged_purity_df['ClusterPurity_db']
    cluster_purity_ms1 = merged_purity_df['ClusterPurity_ms1']

    plt.figure(figsize=(10, 6))

    plt.scatter(cluster_purity_db, cluster_purity_ms1, marker='o', color='blue', label='DF1 vs. DF2')

    # Annotate each point with common cluster indices, avoiding overlaps
    for i, cluster_index in enumerate(common_cluster_indices):
        x = cluster_purity_db.iloc[i]
        y = cluster_purity_ms1.iloc[i]

        # Offset the annotation to avoid overlaps
        offset = 0.01
        x_offset = offset if i % 2 == 0 else -offset
        y_offset = offset if i % 2 == 0 else -offset

        plt.annotate(cluster_index, (x + x_offset, y + y_offset))

    plt.xlabel('Cluster Purity in DF1')
    plt.ylabel('Cluster Purity in DF2')
    plt.title('Scatter Plot of Cluster Purity')
    plt.legend()

    plt.grid(True)
    plt.show()

    x= cluster_purity_db
    y= cluster_purity_ms1

    heatmap, xedges, yedges = np.histogram2d(x, y, bins=20, range=[[0, 1], [0, 1]])

    # Transpose the heatmap
    heatmap = heatmap.T

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the heatmap
    im = ax.imshow(heatmap, extent=[0, 1, 0, 1], cmap='plasma', norm=LogNorm(), origin='lower')

    # Add a colorbar
    cbar = plt.colorbar(im)

    # Set axis labels
    ax.set_xlabel('db results')
    ax.set_ylabel('ms1 results')

    # Show the plot
    plt.show()

    tolerance = 0.05
    # List clusters not within the tolerance
    clusters_not_within_tolerance = []
    for i in range(len(common_cluster_indices)):
        purity1 = cluster_purity_db.iloc[i]
        purity2 = cluster_purity_ms1.iloc[i]
        if abs(purity1 - purity2) > tolerance:
            clusters_not_within_tolerance.append(common_cluster_indices.iloc[i])

    cluster_peptide_portion = merged_data.groupby('#ClusterIdx')['Peptide'].value_counts(normalize=True)

    num_cols = 5

    num_clusters = len(clusters_not_within_tolerance[0:100])
    num_rows = int(np.ceil(num_clusters / num_cols))

    # Create subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20 , 5 * num_rows))

    # Flatten the axes array and iterate through the first 100 clusters to plot the pie chart
    for i, cluster_idx in enumerate(clusters_not_within_tolerance[0:100]):
        portion_data = cluster_peptide_portion.loc[cluster_idx]
        unique_peptide_types_cluster = portion_data.index
        spectra_number = len(merged_data[merged_data['#ClusterIdx'] == cluster_idx]['Peptide'])
        row = i // num_cols
        col = i % num_cols
        ax = axes[row, col]
        ax.pie(portion_data, labels=unique_peptide_types_cluster, autopct='%1.1f%%')
        ax.set_title(f'Cluster {cluster_idx} with {spectra_number} spectra Peptide Fraction')

    # Hide empty subplots if any
    for i in range(num_clusters, num_rows * num_cols):
        fig.delaxes(axes.flatten()[i])

    plt.tight_layout()
    plt.show()

    clusters_within_tolerance = []

    for i in range(len(common_cluster_indices)):
        purity1 = cluster_purity_db.iloc[i]
        purity2 = cluster_purity_ms1.iloc[i]
        if abs(purity1 - purity2) <= tolerance:
            clusters_within_tolerance.append(common_cluster_indices.iloc[i])

    num_clusters = len(clusters_within_tolerance[0:100])
    num_rows = int(np.ceil(num_clusters / num_cols))

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))

    for i, cluster_idx in enumerate(clusters_within_tolerance[0:100]):
        portion_data = cluster_peptide_portion.loc[cluster_idx]
        unique_peptide_types_cluster = portion_data.index
        spectra_number = len(merged_data[merged_data['#ClusterIdx'] == cluster_idx]['Peptide'])
        row = i // num_cols
        col = i % num_cols
        ax = axes[row, col]
        ax.pie(portion_data, labels=unique_peptide_types_cluster, autopct='%1.1f%%')
        ax.set_title(f'Cluster {cluster_idx} with {spectra_number} spectra Peptide Fraction')
        # Hide empty subplots if any
    for i in range(num_clusters, num_rows * num_cols):
        fig.delaxes(axes.flatten()[i])

    plt.tight_layout()
    plt.show()





