import os
from itertools import combinations
from pyteomics import mzml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymzml
from  xml.etree.ElementTree import ParseError

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

def calculate_cluster_purity(cluster, matching_pairs_set):
    total_spectra = len(cluster)
    if total_spectra == 0:
        return 1
    if total_spectra ==1:
        return 1

    max_matching_fraction = 0.0

    for index, spec1 in cluster.iterrows():
        matching_count = 0
        for _, spec2 in cluster.iterrows():
            if spec1['#Filename'] != spec2['#Filename']:
                continue

            if are_matching_spectra(spec1['#Filename'], spec1['#Scan'], spec2['#Scan'], matching_pairs_set):
                matching_count += 1

        matching_fraction = matching_count / total_spectra

        if matching_fraction > max_matching_fraction:
            max_matching_fraction = matching_fraction

    if max_matching_fraction == 0:
        return 1/total_spectra

    return max_matching_fraction




def compare_scans(scan1, scan2, mass_tolerance=0.01, rt_tolerance=60):
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
