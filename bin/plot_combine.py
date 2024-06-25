from matplotlib import pyplot as plt


falcon_x = [96,145,203,273,362,]
falcon_y = [0.9764915441899594,0.9722627486949142,0.9675806275302854,0.961065092182063,0.9531696329137179,]

maracluster_x = [25,32,43,59,82,128]
maracluster_y = [0.9967063573235393,0.995127108862351,0.9928244331292999,0.9897069079966232, 0.9852947966132841,0.9770505998545196]

mscluster_x_MSRT = [32,53,99,148,181,586]
mscluster_y_MSRT = [0.951919117750166,0.9399324174740764,0.9169879561130767,0.8907100412354894,0.8754405110171304,0.6144343047418329]

falcon_x_DB = [96, 145,203,273,362]
falcon_y_DB = [0.9842243711050124,0.9784040298187534,0.972911202042458,0.9662733617461768,0.9580752827454067]


maracluster_x_DB = [25,32,43,59, 82,128]
maracluster_y_DB = [0.9952236308907524,0.993074744759776,0.990518571340414,0.9875523679573239, 0.9837496485947217,0.9778912940627935]

mscluster_x_DB = [32,53,99,148,181,586]
mscluster_y_DB = [0.990518571340414,0.9859990310224227,0.9688700579915228,0.954755816184904,0.9461295812542946,0.8446518259382786]



plt.figure(figsize=(10, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
plt.plot(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
plt.plot(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
plt.plot(falcon_y_DB, falcon_x_DB, '^-', color='red', label='Falcon_DB')
plt.plot(maracluster_y_DB, maracluster_x_DB, '^-', color='green', label='MaraCluster_DB')
plt.plot(mscluster_y_DB, mscluster_x_DB, '^-', color='orange', label='MSCluster_DB')
plt.plot(mscluster_y_MSRT, mscluster_x_MSRT, 'o-', color='orange', label='MSCluster_MSRT')

# Adjust axis labels
plt.xlabel('Average Purity', fontsize=12)
plt.ylabel('N50 Value', fontsize=12)

# Reverse the x-axis
plt.gca().invert_xaxis()


# Add legend, grid, and layout adjustments
plt.legend(loc='lower right', fontsize=12)
plt.title("Clustering Benchmark Results on Joint Proteomics Dataset")
plt.grid(True)
plt.tight_layout()
plt.show()