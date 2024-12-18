from matplotlib import pyplot as plt

falcon_x = [872,1010,1139,1204,1284]
falcon_y = [0.9974634876688682, 0.9973036156255283,0.9971023999263022,0.99694326501403,0.9969418515298243]

maracluster_x = [1,1,2,6,25,137]
maracluster_y = [ 0.9999951110273501,0.9999929966296334,0.9999892783437573,0.9999823528290308,0.9999611395564981,0.9998653171805034]

mscluster_x = [1,1,1,29,284,403,615,720]
mscluster_y = [0.9992571484644659,0.9989708676786451,0.9989504060407933, 0.9985090731730711,0.9979890452122079,0.9977217285260128,0.9976506362730488,0.9975694139661626]







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
plt.legend(loc='upper left', fontsize=14)
plt.title("MSV000081981 Sample Rate: 30%", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()