from matplotlib import pyplot as plt

falcon_x = [322,373,418,437,468]
falcon_y = [0.9926355393861254,0.9921841195396454,0.991933814052273,0.991597514141108,0.991378839262438]

maracluster_x = [1,1,2,5,16,78]
maracluster_y = [0.9999985697425275,0.9999937747826139,0.9999873030564583,0.9999611856268771,0.998999504934814,0.9996744237107452]

mscluster_x = [1,1,1,31,135,181,288,365]
mscluster_y = [0.9989022267684392,0.9991100934985608,0.998269296436166,0.9980544747452545,0.9977748738069953,0.9975114128796164,0.9971564284318974,0.9967185760497436]






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
plt.title("MSV000081981 Sample Rate: 10%", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()