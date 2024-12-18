from matplotlib import pyplot as plt

mscluster_per = [6.548007160443321,6.80620501073829,8.120317821903257,11.181950242852215,15.463775435687166,17.598765476104955,24.21760567512477,28.895803412086877]
falcon_per = [43.469680077782904,49.67984563333054,54.74692132195106,59.25115311616929,63.40533260361594]
maracluster_per= [4.618771102533154,7.575296433981506,11.384742947359994,16.84263421398451,27.791685425997205]

falcon_x = [1867,2223,2355,2680,2786]
falcon_y = [0.9987320374697576,0.9986978636355741,0.9986528010951592,0.9986592417208241,0.9985245444982306]

maracluster_x = [164,25,6,2,1,1]
maracluster_y = [0.9999213432751548,0.9999767562351078,0.9999902112540489,0.99999418455864,0.9999973980404087]

mscluster_x = [1,1,1,13,453,747,1117,1316]
mscluster_y = [0.9995451740539939, 0.9995099382784811,0.9994819886939436,0.9991393755658139,0.9988074937283313,0.9986761814110127,0.9983534011339824,0.9984081490598511]

maracluster_y =maracluster_y[::-1]







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
plt.legend(loc='upper right', fontsize=14)
plt.title("MSV000081981 Sample Rate: 70%", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()