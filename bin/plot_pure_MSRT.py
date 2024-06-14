# Complete code to generate the requested plot with purity as line plots and count as bar plots

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Data for each method, arranged in the order of bins
data = {
    "Bin": ['2 to 2', '3 to 4', '5 to 8', '9 to 16', '17 to 32',
            '33 to 64', '65 to 128', '129 to 256', '257 to 512',
            '513 to 1024', '1025 to 2048', '2049 to 4096', '4097 to 8192', '8193 to 16384'],
    "mscluster_Purity": [0.996015, 0.991748, 0.984973, 0.978324, 0.979238,
                         0.977536, 0.975906, 0.975667, 0.961718, 0.958553,
                         0.940228, 0.934496, 0.936255, 0.589451],
    "mscluster_Count": [16562, 8907, 5328, 2540, 2024,
                        1662, 1233, 837, 485, 251, 147, 139, 55, 2],
    "falcon_Purity": [0.999084, 0.997562, 0.994193, 0.989909, 0.986410,
                      0.982499, 0.980914, 0.984315, 0.979776, 0.969567,
                      0.971448, 0.961412, 0.934286, 0.787866],
    "falcon_Count": [57318, 26215, 11188, 4977, 2701,
                     1696, 1203, 772, 579, 446, 339, 219, 76, 6],
    "maracluster_Purity": [0.999414, 0.999013, 0.998330, 0.998130,
                           0.998085, 0.997733, 0.998906, 0.998888, 0.998147, 0.999107,
                           0.999751, 0.996158, 1.000000,0],
    "maracluster_Count": [88671, 50420, 25568, 12505,
                          5990, 3055, 1356, 664, 318, 156, 60, 22, 1,0]
}

# Create a DataFrame
df = pd.DataFrame(data)

# Reordering based on the provided bin order
df = df.set_index('Bin')
df = df.reindex(data['Bin'])
# Plot adjustments
plt.rcParams['font.size'] = 16  # Adjusts the global font size
plt.rcParams['legend.fontsize'] = 10  # Adjusts the legend font size
# Start plotting
fig, ax1 = plt.subplots(figsize=(15, 7))

# Set positions for the bars
positions = np.arange(len(df.index))

# Bar plots for Cluster Count
bar_width = 0.2  # width of the bars
ax1.bar(positions - bar_width, df['mscluster_Count'], width=bar_width, label='mscluster Count')
ax1.bar(positions, df['falcon_Count'], width=bar_width, label='falcon Count')
ax1.bar(positions + bar_width, df['maracluster_Count'], width=bar_width, label='maracluster Count')

# Set labels and title
ax1.set_xlabel('Cluster Size Range')
ax1.set_ylabel('Cluster Count (log scale)')
ax1.set_yscale('log')
ax1.set_title('Average Cluster Purity and Count by Size Range')
ax1.set_xticks(positions)
ax1.set_xticklabels(df.index, rotation=45)
ax1.tick_params(axis='x', labelsize=14)  # Increase X tick label size
ax1.tick_params(axis='y', labelsize=14)  # Increase Y tick label size for ax1
ax1.legend(loc='upper left')

# Create another y-axis for the Purity
ax2 = ax1.twinx()
ax2.plot(positions, df['mscluster_Purity'], marker='o', label='mscluster Purity')
ax2.plot(positions, df['falcon_Purity'], marker='o', label='falcon Purity')
ax2.plot(positions, df['maracluster_Purity'], marker='o', label='maracluster Purity')
ax2.tick_params(axis='y', labelsize=14)
ax2.set_ylabel('Average Purity')
ax2.legend(loc='lower left')

# Show grid and plot
ax1.grid(True)
plt.tight_layout()
plt.show()
