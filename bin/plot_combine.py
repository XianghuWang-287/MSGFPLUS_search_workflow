from matplotlib import pyplot as plt


falcon_x = [8, 24, 25, 26, 26]
falcon_y = [0.9896683397928712, 0.9748935149079129, 0.973328106397189, 0.972145244693312, 0.9704464146834656]

# maracluster_x = [17,17,16,14,11]
maracluster_x = [25, 23, 21, 18, 16, 14]
maracluster_y = [0.972550831792976, 0.9771322946923686, 0.9806311064166887, 0.9839054660681278, 0.9875231053604436,
                 0.9911407446527595]

falcon_x_DB = [8, 24, 25, 26, 26]
falcon_y_DB = [0.995044158, 0.98821705, 0.986540639, 0.984972717, 0.983399087]

# maracluster_x = [17,17,16,14,11]
maracluster_x_DB = [25, 23, 21, 18, 16, 14]
maracluster_y_DB = [0.9845260100343279, 0.98844732, 0.990071297, 0.991444415, 0.992632691, 0.993992606]

mscluster_x_DB = [17, 21, 25, 29, 30, 37, 40]
mscluster_y_DB = [0.9908740620563781, 0.9879548643243082, 0.9824862146795806, 0.9740523752382297, 0.9682559792784042,
                  0.9473964314725077, 0.9362784845519363]

mscluster_x_MSRT = [17, 21, 25, 29, 30, 37, 40]
mscluster_y_MSRT = [0.9776076522679645, 0.9692229618794591, 0.9601011399674295, 0.9533140396696548, 0.9497304221743599,
                    0.9363323156214112, 0.9283237802320704]


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