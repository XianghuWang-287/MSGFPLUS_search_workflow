from matplotlib import pyplot as plt

mscluster_per = [6.1730665035690375,6.6252878788072245,8.155392202554317,12.003849308346028,16.688862755813652,19.007300661881846]
falcon_per = [34.02837098581706,40.420759536125146,45.534670376371615,50.06497091099305,54.1512645386512]
maracluster_per= [4.5134885399655165,7.586636039042985,11.521791962845066,17.070090695390917,27.697229170470646]


falcon_x = [872,1010,1139,1204,1284]
falcon_y = [0.9974634876688682, 0.9973036156255283,0.9971023999263022,0.99694326501403,0.9969418515298243]

maracluster_x = [1,1,2,6,25,137]
maracluster_y = [0.9999929966296334,0.9999892783437573,0.9999823528290308,0.9999611395564981,0.9998653171805034]

mscluster_x = [1,1,1,29,284,403,615,720]
mscluster_y = [0.9989504060407933, 0.9985090731730711,0.9979890452122079,0.9977217285260128,0.9976506362730488,0.9975694139661626]










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
plt.title("MSV000081981 Sample Rate: 30%", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()