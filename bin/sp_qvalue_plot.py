import pandas as pd
import matplotlib.pyplot as plt

# Load the TSV file into a DataFrame
df = pd.read_csv('../data/results/output.tsv')
print("Column names:", df.columns)

# Extract the qvalue column
qvalues = df['QValue']
# Initialize a dictionary to store the count of spectra for each qvalue
qvalue_counts = {qval: sum(qvalues <= qval) for qval in qvalues}

# Sort the dictionary by qvalue
sorted_qvalue_counts = dict(sorted(qvalue_counts.items()))

# Extract the qvalues and counts for plotting
qvalue_thresholds = list(sorted_qvalue_counts.keys())
spectra_counts = list(sorted_qvalue_counts.values())

# Plot the figure
plt.plot(qvalue_thresholds, spectra_counts, marker='o')
plt.xlabel('qvalue')
plt.ylabel('Number of Spectra')
plt.title('Number of Spectra vs. qvalue')
plt.grid(True)
plt.show()
