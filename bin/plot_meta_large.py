from matplotlib import pyplot as plt


mscluster_x = [1,10,471,1743]
mscluster_y = [0.9987760692342341,0.9977736639071505,0.996543148551809,0.9915985003865873]

falcon_x = [2877,4785,5127]
falcon_y = [0.99849630154489,0.9981019744016099,0.9980665816962624]

maracluster_x = [1,2,6,185]
maracluster_y = [0.9999958737481095,0.9999868134694759,0.9999930406241849,0.9999047623371281]



plt.figure(figsize=(8, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
plt.plot(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
plt.plot(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
plt.plot(mscluster_y, mscluster_x, 'o-', color='orange', label='MSCluster_MSRT')

# Adjust axis labels
plt.xlabel('Average Purity', fontsize=14)
plt.ylabel('N10 Value', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# Reverse the x-axis
plt.gca().invert_xaxis()
plt.xlim(right = 0.99)


# Add legend, grid, and layout adjustments
plt.legend(loc='lower right', fontsize=14)
plt.title("Clustering Benchmark Results on Metabolomics Dataset",fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()