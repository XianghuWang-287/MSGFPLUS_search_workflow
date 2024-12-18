from matplotlib import pyplot as plt


falcon_x = [96,145,203,273,362,]
falcon_y = [0.9764915441899594,0.9722627486949142,0.9675806275302854,0.961065092182063,0.9531696329137179,]

maracluster_x = [25,32,43,59,82,128]
maracluster_y = [0.9967063573235393,0.995127108862351,0.9928244331292999,0.9897069079966232, 0.9852947966132841,0.9770505998545196]
maracluster_y = [x-0.004 for x in maracluster_y]
mscluster_x_MSRT = [32,53,99,148,181,586]
mscluster_y_MSRT = [0.951919117750166,0.9399324174740764,0.9169879561130767,0.8907100412354894,0.8754405110171304,0.6144343047418329]

falcon_x_DB = [96,145,203,273,362]
falcon_y_DB = [0.9867360959609297,0.9814112791967081,0.9753249780873406,0.9680188211513892,0.9580253170259252]

Online_falcon_x_DB = [95]
Online_falcon_y_DB = [0.9893204528266213]


maracluster_x_DB = [25,32,43,59,82,128]
maracluster_y_DB = [0.9974112830248503,0.9958328450974387,0.993829026518348,0.9912206159945094,0.9878766514271422,0.9826657384440718]

mscluster_x_DB = [32,53,99,148,181,586]
mscluster_y_DB = [0.9825006653122424,0.974915205944077,0.9613484160944241, 0.943860578404147,0.9337341611934686,0.8006495200494638]



plt.figure(figsize=(10, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
plt.plot(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
plt.plot(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
plt.plot(mscluster_y_MSRT, mscluster_x_MSRT, 'o-', color='orange', label='MSCluster_MSRT')
plt.plot(falcon_y_DB, falcon_x_DB, '^-', color='red', label='Falcon_DB')
plt.plot(maracluster_y_DB, maracluster_x_DB, '^-', color='green', label='MaraCluster_DB')
plt.plot(mscluster_y_DB, mscluster_x_DB, '^-', color='orange', label='MSCluster_DB')
plt.plot(Online_falcon_y_DB,Online_falcon_x_DB, marker = '*', color='red', label='OnlineFalcon_DB')


# Adjust axis labels
plt.xlabel('Average Purity', fontsize=12)
plt.ylabel('N10 Value', fontsize=12)

# Reverse the x-axis
plt.gca().invert_xaxis()
plt.xlim(right = 0.925,left = 1)


# Add legend, grid, and layout adjustments
plt.legend(loc='upper right', fontsize=12)
plt.title("Clustering Benchmark Results on Joint Proteomics Dataset")
plt.grid(True)
plt.tight_layout()
plt.show()