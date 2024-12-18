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
from collections import defaultdict
from sklearn.metrics import r2_score

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
        node_a =f"{spec1['#Filename']}_{spec1['#Scan']}"
        G.add_node(node_a, filename=spec1['#Filename'])
        for _, spec2 in cluster.iterrows():
            if spec1['#Filename'] != spec2['#Filename']:
                continue
            if are_matching_spectra(spec1['#Filename'], spec1['#Scan'], spec2['#Scan'], matching_pairs_set):
                node_a =f"{spec1['#Filename']}_{spec1['#Scan']}"
                node_b =f"{spec2['#Filename']}_{spec2['#Scan']}"
                G.add_node(node_a,filename=spec1['#Filename'])
                G.add_node(node_b,filename=spec2['#Filename'])
                G.add_edge(node_a,node_b)

    return G

def create_matching_network_across_files(cluster,tolerence):
    G = nx.Graph()

    for _, spec1 in cluster.iterrows():
        node_a =f"{spec1['#Filename']}_{spec1['#Scan']}"
        G.add_node(node_a, filename=spec1['#Filename'])
        for _, spec2 in cluster.iterrows():
            if (spec1['#Filename'] == spec2['#Filename'] and spec1['#Scan'] == spec2['#Scan']):
                continue
            if abs(spec1['']):
                node_a =f"{spec1['#Filename']}_{spec1['#Scan']}"
                node_b =f"{spec2['#Filename']}_{spec2['#Scan']}"
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


    for node in G.nodes:
        if node not in [c for component in components for c in component]:
            components.append({node})

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster['#Filename'])

    max_fraction = 0.0

    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max((comp.number_of_nodes() for comp in S if any(comp.nodes[node]['filename'] == filename for node in comp.nodes())), default=0)
        fraction = largest_component_size / count
        max_fraction = max(max_fraction, fraction)

    return max_fraction

def calculate_cluster_purity_weighted_avg(cluster, matching_pairs_set):
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, matching_pairs_set)

    # Get connected components
    components = list(nx.connected_components(G))

    for node in G.nodes:
        if node not in [c for component in components for c in component]:
            components.append({node})

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster['#Filename'])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max((comp.number_of_nodes() for comp in S if any(comp.nodes[node]['filename'] == filename for node in comp.nodes())), default=0)
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    weighted_sum = sum(value * frequency for value, frequency in zip(values, frequencies))

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    weighted_average = weighted_sum / total_frequency
    return weighted_average

def calculate_cluster_purity_avg(cluster, matching_pairs_set):
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, matching_pairs_set)

    # Get connected components
    components = list(nx.connected_components(G))

    for node in G.nodes:
        if node not in [c for component in components for c in component]:
            components.append({node})

    S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    # Calculate the count of each filename in the matching pairs within the cluster
    file_counts = Counter(cluster['#Filename'])

    frequencies = []
    values = []
    # Calculate the fraction of the largest component for each file
    for filename, count in file_counts.items():
        largest_component_size = max((comp.number_of_nodes() for comp in S if any(comp.nodes[node]['filename'] == filename for node in comp.nodes())), default=0)
        fraction = largest_component_size / count
        frequencies.append(count)
        values.append(fraction)
    total_purity = sum(values)

    # Calculate the total frequency
    total_frequency = sum(frequencies)

    # Calculate the weighted average
    average = total_purity / total_frequency
    return average

def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=30):
    mass_diff = abs(scan1[0] - scan2[0])
    rt_diff = abs(scan1[1] - scan2[1])
    return mass_diff <= mass_tolerance and rt_diff <= rt_tolerance

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
           and abs(spec1[1] - spec2[1]) <= 0.01 and abs(spec1[2] - spec2[2]) <= 1
    ]
    G.add_edges_from(edges)

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
def calculate_cluster_purity_weighted_avg_new(cluster, method):
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

if __name__ == "__main__":
    #constructing the ms1 results
    database_results = pd.read_csv('./PXD023047_filtered.tsv', sep='\t')  # Adjust file path and format accordingly
    folder_path = '/home/user/LabData/XianghuData/MS_Cluster_datasets/PXD023047_convert/mzML'

    #cluster_results = pd.read_csv('/home/user/LAB_share/XianghuData/MS_Cluster_datasets/PXD023047_results/msculster_results/clustering/mscluster_clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly

    cluster_results = pd.read_csv('../data/results/nf_output/clustering/clusterinfo.tsv',sep='\t')  # Adjust file path and format accordingly

    cluster_sizes = cluster_results.groupby('#ClusterIdx').size()

    clusters_filter_index = cluster_sizes[cluster_sizes >= 7].index

    cluster_results = cluster_results[cluster_results['#ClusterIdx'].isin(clusters_filter_index)]
    #handle the new version workflow filename issue
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzML')
    # matching_pairs_all_files = []
    #
    #
    # for filename in os.listdir(folder_path):
    #     if filename.endswith('.mzML'):
    #         mzml_file = os.path.join(folder_path, filename)
    #         ms1_data = get_ms1_data(mzml_file)
    #
    #         matching_pairs = []
    #
    #         for scan1, scan2 in combinations(ms1_data, 2):
    #             if compare_scans(scan1, scan2):
    #                 matching_pairs.append((scan1[2], scan2[2]))
    #
    #         matching_pairs_all_files.extend([(filename, pair[0], pair[1]) for pair in matching_pairs])
    #
    # matching_pairs_all_files = [(f'mzML/{filename}', scan_start, scan_end) for filename, scan_start, scan_end in matching_pairs_all_files]
    #
    # matching_pairs_set = {(item[0], item[1], item[2]) for item in matching_pairs_all_files}

    cluster_purity = cluster_results.groupby('#ClusterIdx').apply(lambda x: calculate_cluster_purity_weighted_avg_new(x,'mscluster'))

    cluster_indices = cluster_results['#ClusterIdx'].unique()

    purity_by_cluster_index = pd.DataFrame({'Cluster_indices': cluster_indices, 'ClusterPurity': cluster_purity})

    #constructing the db results
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
    # Filter database results based on "DB:EValue" column
    #filtered_database_results = database_results[database_results['DB:EValue'] < 0.002]
    filtered_database_results = database_results

    # Extract unique spectra
    unique_spectra_db = filtered_database_results[['MzIDFileName', 'ScanNumber']].drop_duplicates()

    # Filter cluster_results to only include spectra present in filtered database results
    filtered_cluster_results = cluster_results.merge(unique_spectra_db,
                                                     left_on=['#Filename', '#Scan'],
                                                     right_on=['MzIDFileName', 'ScanNumber'],
                                                     how='inner')
    # Merge the filtered clusters with the corresponding peptide identifications based on filename, scan number, and charge
    merged_data = pd.merge(filtered_cluster_results, database_results,
                           left_on=['#Filename', '#Scan', '#Charge'],
                           right_on=['MzIDFileName', 'ScanNumber', 'Charge'],
                           how='inner')

    db_sizes = merged_data.groupby('#ClusterIdx').size()

    merged_filter_index = db_sizes[db_sizes >= 7].index
    # merged_filter_index = db_sizes.index

    merged_data = merged_data[merged_data['#ClusterIdx'].isin(merged_filter_index)]
    #handle the "I" and "L" replacement case
    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')

    # cluster_purity = merged_data.groupby('#ClusterIdx').apply(lambda x: calculate_cluster_purity_weighted_avg(x, matching_pairs_set))
    #
    # cluster_indices = merged_data['#ClusterIdx'].unique()
    #
    # purity_by_cluster_index = pd.DataFrame({'Cluster_indices': cluster_indices, 'ClusterPurity': cluster_purity})



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

    cluster_sizes = merged_data.groupby('#ClusterIdx').size().loc[common_cluster_indices]
    weights = cluster_sizes*2

    plt.figure(figsize=(10, 6))

    plt.scatter(cluster_purity_db, cluster_purity_ms1, marker='o', color='blue', label='db results vs. ms1 results')

    # Annotate each point with common cluster indices, avoiding overlaps
    for i, cluster_index in enumerate(common_cluster_indices):
        x = cluster_purity_db.iloc[i]
        y = cluster_purity_ms1.iloc[i]

        # Offset the annotation to avoid overlaps
        offset = 0.01
        x_offset = offset if i % 2 == 0 else -offset
        y_offset = offset if i % 2 == 0 else -offset

        plt.annotate(cluster_index, (x + x_offset, y + y_offset))

    plt.xlabel('Cluster Purity in db results')
    plt.ylabel('Cluster Purity in ms1 results')
    plt.title('Scatter Plot of Cluster Purity')
    plt.legend()

    plt.grid(True)
    plt.show()

    x= cluster_purity_db
    y= cluster_purity_ms1
    y_lower = x - 0.1
    y_upper = x + 0.1
    y_pred = np.clip(y, y_lower, y_upper)  # Clip values to ensure they are within ±0.1 range
    r2 = r2_score(y, y_pred)
    print(f'R^2 between MS1 and DB purity (y=x): {r2:.4f}')


    heatmap, xedges, yedges = np.histogram2d(x, y, bins=20, range=[[0, 1], [0, 1]])

    # Transpose the heatmap
    heatmap = heatmap.T

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the heatmap


    # Now use LogNorm with the offset data
    # vmin should be the smallest value you've added as an offset, to ensure 0s are included and distinctly represented
    im = ax.imshow(heatmap, extent=[0, 1, 0, 1], cmap='plasma', origin='lower')

    # Add a colorbar
    cbar = plt.colorbar(im)

    # Set axis labels
    ax.set_xlabel('db results')
    ax.set_ylabel('ms1 results')
    plt.tight_layout()
    # Show the plot
    plt.show()

    tolerance = 0.1
    # List clusters not within the tolerance
    clusters_not_within_tolerance = []
    for i in range(len(common_cluster_indices)):
        purity1 = cluster_purity_db.iloc[i]
        purity2 = cluster_purity_ms1.iloc[i]
        if abs(purity1 - purity2) > tolerance:
            clusters_not_within_tolerance.append(common_cluster_indices.iloc[i])

    below_y_equals_x_minus_tolerance = 0
    above_y_equals_x_plus_tolerance = 0

    for i in range(len(common_cluster_indices)):
        purity_db = cluster_purity_db.iloc[i]
        purity_ms1 = cluster_purity_ms1.iloc[i]

        # Check if the point is within tolerance
        if abs(purity_db - purity_ms1) > tolerance:
            clusters_not_within_tolerance.append(common_cluster_indices.iloc[i])

        # Check if the point is below y = x - 0.1
        if purity_ms1 < purity_db - tolerance:
            below_y_equals_x_minus_tolerance += 1

        # Check if the point is above y = x + 0.1
        if purity_ms1 > purity_db + tolerance:
            above_y_equals_x_plus_tolerance += 1

    print("Clusters not within tolerance:", clusters_not_within_tolerance)
    print("Points below y = x - 0.1:", below_y_equals_x_minus_tolerance)
    print("Points above y = x + 0.1:", above_y_equals_x_plus_tolerance)

    cluster_peptide_portion = merged_data.groupby('#ClusterIdx')['Peptide'].value_counts(normalize=True)

    num_cols = 5

    num_clusters = len(clusters_not_within_tolerance[0:100])
    num_rows = int(np.ceil(num_clusters / num_cols))

    # Create subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(25 , 5 * num_rows))

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
    print(len(clusters_within_tolerance))
    print(len(clusters_not_within_tolerance))





