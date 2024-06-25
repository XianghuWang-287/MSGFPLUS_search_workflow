from matplotlib import pyplot as plt


falcon_x = [299,448,512,522,557]
falcon_y = [0.9333594902993452,0.9282725146345037,0.9240639595317547,0.921489199831977,0.9175995705688137]

maracluster_x = [1,1,1,1,6]
maracluster_y = [0.9996361978668292,0.9988184890445627,0.9968141221685878,0.9948826748809761,0.9886893665160859]

mscluster_x_MSRT = [1,1,1,4,6,23,40]
mscluster_y_MSRT = [0.9495690763254434,0.9442609560901015,0.9435251501999243,0.9373086296542512,0.9272055662910907,0.9004527723056273,0.89207518585509]



plt.figure(figsize=(10, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
plt.plot(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
plt.plot(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
plt.plot(mscluster_y_MSRT, mscluster_x_MSRT, 'o-', color='orange', label='MSCluster_MSRT')

# Adjust axis labels
plt.xlabel('Average Purity', fontsize=12)
plt.ylabel('N50 Value', fontsize=12)

# Reverse the x-axis
plt.gca().invert_xaxis()


# Add legend, grid, and layout adjustments
plt.legend(loc='upper left', fontsize=12)
plt.title("Clustering Benchmark Results on Metabolomics Dataset")
plt.grid(True)
plt.tight_layout()
plt.show()