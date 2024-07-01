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
import copy

method_dic = {'mscluster':{'filename':'#Filename','scan':'#Scan','mass':'#ParentMass','rt_time':'#RetTime'},
              'falcon':{'filename':'filename','scan':'scan','mass':'precursor_mz','rt_time':'retention_time'},
              'maracluster':{'filename':'filename','scan':'scan','mass':'#ParentMass','rt_time':'#RetTime'}}

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
            if abs(spec1[method_dic[method]['mass']]-spec2[method_dic[method]['mass']])<=0.01 and abs(spec1[method_dic[method]['rt_time']]-spec2[method_dic[method]['rt_time']])<=0.5:
                node_a = f"{spec1[method_dic[method]['filename']]}_{spec1[method_dic[method]['scan']]}"
                node_b = f"{spec2[method_dic[method]['filename']]}_{spec2[method_dic[method]['scan']]}"
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

def calculate_cluster_purity_weighted_avg(cluster,method):
    if len(cluster) == 1:
        return 1
    # Create a matching network considering only matching pairs with the same filename and scan number
    G = create_matching_network(cluster, method)

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

def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=65):
    mass_diff = abs(scan1[0] - scan2[0])
    rt_diff = abs(scan1[1] - scan2[1])
    return mass_diff <= mass_tolerance and rt_diff <= rt_tolerance

def mscluster_merge(cluster_results,database_results):
    cluster_results['#Filename'] = cluster_results['#Filename'].str.replace('input_spectra', 'mzML')
    cluster_results['#RetTime'] = cluster_results['#RetTime']/60
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
    reference_data['#RetTime'] = reference_data['#RetTime'] / 60
    database_results['MzIDFileName'] = 'mzML/' + database_results['MzIDFileName']
    merged_data = pd.merge(cluster_results, reference_data, left_on=['filename', 'scan'], right_on=['#Filename', '#Scan'],how='inner')
    merged_data = merged_data.merge(database_results, left_on=['filename', 'scan'],
                                        right_on=['MzIDFileName', 'ScanNumber'], how='inner')
    return merged_data
if __name__ == "__main__":
    #constructing the ms1 results
    database_results = pd.read_csv('./Combine_test_filtered.tsv', sep='\t')  # Adjust file path and format accordingly
    print("number of identified peptides:", len(database_results))
    database_results = database_results[database_results['DB:EValue'] < 0.002]

    mscluster_results =  pd.read_csv('../data/Combine_results/mscluster_clusterinfo.tsv',sep='\t')
    falcon_results = pd.read_csv('/home/user/LabData/XianghuData/Falcon_Cluster_Benchmark/combine_0.3/output_summary/cluster_info.tsv', sep='\t')
    maracluster_results = pd.read_csv('../data/Combine_results/MaRaCluster_processed.clusters_p10.tsv', sep='\t')


    cluster_results =  falcon_merge(falcon_results,database_results)
    cluster_sizes = cluster_results.groupby('cluster').size()
    clusters_filter_index = cluster_sizes[cluster_sizes >= 10].index
    cluster_results = cluster_results[cluster_results['cluster'].isin(clusters_filter_index)]
    cluster_purity = cluster_results.groupby('cluster').apply(lambda x: calculate_cluster_purity_weighted_avg(x, 'falcon'))
    print("finished calculating cluster_purity")

    cluster_indices = cluster_results['cluster'].unique()

    purity_by_cluster_index = pd.DataFrame({'Cluster_indices': cluster_indices, 'ClusterPurity': cluster_purity})

    merged_data = cluster_results

    merged_data['Peptide'] = merged_data['Peptide'].str.replace('I', 'L')

    merged_data.to_csv("../data/Combine_results/merged_data.csv", index=False)

    cluster_purity_db = merged_data.groupby('cluster')['Peptide'].apply(
        lambda x: x.value_counts().max() / len(x))

    # Calculate the number of spectra in each cluster
    cluster_indices_db = merged_data['cluster'].unique()

    purity_by_cluster_index_db = pd.DataFrame({'Cluster_indices': cluster_indices_db, 'ClusterPurity': cluster_purity_db})

    merged_purity_df = purity_by_cluster_index_db.merge(purity_by_cluster_index, on='Cluster_indices', suffixes=('_db', '_ms1'), how='inner')

    # Extract the common cluster indices and ClusterPurity values
    common_cluster_indices = merged_purity_df['Cluster_indices']
    cluster_purity_db = merged_purity_df['ClusterPurity_db']
    cluster_purity_ms1 = merged_purity_df['ClusterPurity_ms1']

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

    heatmap, xedges, yedges = np.histogram2d(x, y, bins=20, range=[[0, 1], [0, 1]])

    # Transpose the heatmap
    heatmap = heatmap.T

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the heatmap
    im = ax.imshow(heatmap, extent=[0, 1, 0, 1], cmap='plasma_r', norm=LogNorm(), origin='lower')

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

    cluster_peptide_portion = merged_data.groupby('cluster')['Peptide'].value_counts(normalize=True)

    num_cols = 5

    num_clusters = len(clusters_not_within_tolerance[0:100])
    num_rows = int(np.ceil(num_clusters / num_cols))

    # Create subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20 , 5 * num_rows))

    # Flatten the axes array and iterate through the first 100 clusters to plot the pie chart
    for i, cluster_idx in enumerate(clusters_not_within_tolerance[0:100]):
        portion_data = cluster_peptide_portion.loc[cluster_idx]
        unique_peptide_types_cluster = portion_data.index
        spectra_number = len(merged_data[merged_data['cluster'] == cluster_idx]['Peptide'])
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
        spectra_number = len(merged_data[merged_data['cluster'] == cluster_idx]['Peptide'])
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

    filtered_clusters = merged_purity_df[
        (merged_purity_df['ClusterPurity_ms1'] == 1) & (merged_purity_df['ClusterPurity_db'] != 1)]

    # Extract relevant data for plotting
    cluster_indices_filtered = filtered_clusters['Cluster_indices']
    cluster_purity_db_filtered = filtered_clusters['ClusterPurity_db']
    cluster_purity_ms1_filtered = filtered_clusters['ClusterPurity_ms1']

    # Plot the filtered clusters
    plt.figure(figsize=(10, 6))
    plt.scatter(cluster_purity_db_filtered, cluster_purity_ms1_filtered, marker='o', color='red',
                label='Filtered Clusters')
    plt.xlabel('Cluster Purity in db results')
    plt.ylabel('Cluster Purity in ms1 results')
    plt.title('Filtered Clusters where MS1 Purity is 1 but DB Purity is not')
    plt.legend()
    plt.grid(True)
    plt.show()

    print("Cluster ID numbers where MS1 purity is 1 but DB purity is not 1:")
    for cluster_id in filtered_clusters['Cluster_indices']:
        print(cluster_id)





