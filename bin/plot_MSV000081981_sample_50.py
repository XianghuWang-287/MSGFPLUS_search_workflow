from matplotlib import pyplot as plt

falcon_x = [1367,1582,1762,1927,1984]
falcon_y = [0.9983848155146413,0.9984178518977926,0.9984116537574407,0.9983031543418371,0.9982010295794757]

maracluster_x = [167,27,6,2,1,1]
maracluster_y = [0.9999130082549189,0.9999735122387372,0.9999897225033207,0.9999944035600657,0.9999971756721113,0.9999984044747359]

mscluster_x = [1,1,1,19,421,572]
mscluster_y = [0.9993456220976422,0.9992852855316049,0.9991655795933858,0.9988258110557806,0.9985219186863091,0.9983594629903335]




plt.figure(figsize=(10, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
# plt.scatter(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
# plt.scatter(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
# plt.scatter(mscluster_y, mscluster_x, 'o-', color='orange', label='MSCluster_MSRT')

plt.scatter(falcon_y, falcon_x, color='red', label='Falcon_MSRT')
plt.scatter(maracluster_y, maracluster_x, color='green', label='MaraCluster_MSRT')
plt.scatter(mscluster_y, mscluster_x, color='orange', label='MSCluster_MSRT')

# Adjust axis labels
plt.xlabel('Average Purity', fontsize=14)
plt.ylabel('N10 Value', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
# Reverse the x-axis
plt.gca().invert_xaxis()
plt.xlim(right=0.991)


# Add legend, grid, and layout adjustments
plt.legend(loc='upper right', fontsize=14)
plt.title("MSV000081981 Sample Rate: 50%", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()