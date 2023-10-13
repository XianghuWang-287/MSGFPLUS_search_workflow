import pandas as pd

# Read the CSV file into a DataFrame
df = pd.read_csv('../data/results/MSGF_PLUS-465acd67-group_by_spectrum_old-main.tsv',sep='\t')

# Sort the DataFrame by qvalue in descending order within each group
sorted_df = df.sort_values(by='QValue', ascending=False).groupby('ScanNum').first().reset_index()

# Save the merged DataFrame to a new CSV file
sorted_df.to_csv('../data/results/output.tsv', index=False)

print('Merging completed. Merged data with highest qvalue saved to output.csv.')
