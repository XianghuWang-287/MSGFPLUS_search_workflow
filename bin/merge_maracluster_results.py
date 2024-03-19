import pandas as pd

def process_maracluster_results(input_file, output_file):
    # Read the input file, considering the format described
    # Each spectrum is on a new line, clusters are separated by an empty line
    # Columns are separated by tabs
    # Use on_bad_lines='skip' to skip any malformed lines
    df = pd.read_csv(input_file, sep='\t', header=None, names=['filename', 'scan', 'cluster'], on_bad_lines='skip')

    # Write the processed DataFrame to a new TSV file
    df.to_csv(output_file, sep='\t', index=False)

# Specify the input and output file paths
input_file = '../data/Combine_results/MaRaCluster.clusters_p10.tsv'
output_file = '../data/Combine_results/MaRaCluster_processed.clusters_p10.tsv'

# Call the function with the specified file paths
process_maracluster_results(input_file, output_file)