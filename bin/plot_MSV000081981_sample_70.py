from matplotlib import pyplot as plt

falcon_x = [1867,2223,2355,2680,2786]
falcon_y = [0.9987320374697576,0.9986978636355741,0.9986528010951592,0.9986592417208241,0.9985245444982306]

maracluster_x = [164,25,6,2,1,1]
maracluster_y = [0.9999213432751548,0.9999767562351078,0.9999902112540489,0.99999418455864,0.9999973980404087,0.9999982982454304]

mscluster_x = [1,1,1,13,453,747,1117,1316]
mscluster_y = [0.9995451740539939, 0.9995099382784811,0.9994819886939436,0.9991393755658139,0.9988074937283313,0.9986761814110127,0.9983534011339824,0.9984081490598511]




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
plt.title("MSV000081981 Sample Rate: 70%", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()