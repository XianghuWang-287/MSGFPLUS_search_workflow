from matplotlib import pyplot as plt

mscluster_per = [5.384711696196322,6.13479218516039,8.034331737181812,12.961789461311637,18.08957527379024,21.111744234179447,31.807705476826843,38.80883998793943]
falcon_per = [34.145334961850786,40.83422339648097,46.10203445438703,50.654899094946316,54.73875581947987]
maracluster_per= [3.867175322850178,6.744847428696998,10.49320570935052,15.70339688980422,22.749503012586942]



falcon_y = [0.9926355393861254,0.9921841195396454,0.991933814052273,0.991597514141108,0.991378839262438]


maracluster_y = [0.9999937747826139,0.9999873030564583,0.9999611856268771,0.998999504934814,0.9996744237107452]




mscluster_y = [0.9989022267684392,0.9991100934985608,0.998269296436166,0.9980544747452545,0.9977748738069953,0.9975114128796164,0.9971564284318974,0.9967185760497436]






plt.figure(figsize=(10, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
# plt.scatter(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
# plt.scatter(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
# plt.scatter(mscluster_y, mscluster_x, 'o-', color='orange', label='MSCluster_MSRT')

plt.scatter(falcon_y, falcon_per, color='red', label='Falcon_MSRT')
plt.scatter(maracluster_y, maracluster_per, color='green', label='MaraCluster_MSRT')
plt.scatter(mscluster_y, mscluster_per, color='orange', label='MSCluster_MSRT')

# Adjust axis labels
plt.xlabel('Average Purity', fontsize=14)
plt.ylabel('Percentage of Clustered Spectra', fontsize=14)
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