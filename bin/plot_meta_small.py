from matplotlib import pyplot as plt


falcon_x = [299,448,512,522,557]
falcon_y = [0.9333594902993452,0.9282725146345037,0.9240639595317547,0.921489199831977,0.9175995705688137]

maracluster_x = [1,1,1,1,6]
maracluster_y = [0.9996361978668292,0.9988184890445627,0.9968141221685878,0.9948826748809761,0.9886893665160859]

mscluster_x_MSRT = [1,1,1,4,6,23,40]
mscluster_y_MSRT = [0.9495690763254434,0.9442609560901015,0.9435251501999243,0.9373086296542512,0.9272055662910907,0.9004527723056273,0.89207518585509]

online_falcon_x=[1,1,1,1,1]
online_falcon_y=[0.9430195599075247,0.9379409544974491,0.9369734258740184,0.932700957275679,0.9305002053100746]

online_falcon_x_new=[124,146,152,154,166,167]
online_falcon_y_new=[0.9558041479112528,0.9523845843543103,0.950283718111141,0.9479304669128918,0.9461500292473213,0.938031681235883]

online_falcon_x_sample = [215,224,243,249,259]
online_falcon_y_sample = [0.9385237438642566,0.9351909240338089,0.9317303295380078,0.9304386961924086,0.9290068116213334]


plt.figure(figsize=(8, 6))

# Plot scatter points
# plt.scatter(valid_purities, valid_n50, color='blue', label='Methods_Default_MSRT')
plt.plot(falcon_y, falcon_x, 'o-', color='red', label='Falcon_MSRT')
plt.plot(maracluster_y, maracluster_x, 'o-', color='green', label='MaraCluster_MSRT')
plt.plot(mscluster_y_MSRT, mscluster_x_MSRT, 'o-', color='orange', label='MSCluster_MSRT')
# plt.plot(online_falcon_y, online_falcon_x, 'o-', label='online_falcon_MSRT')
plt.plot(online_falcon_y_new, online_falcon_x_new, 'o-', label='online_falcon_MSRT')
plt.plot(online_falcon_y_sample, online_falcon_x_sample, 'o-', label='online_falcon_sample_MSRT')

# Adjust axis labels
plt.xlabel('Average Purity', fontsize=14)
plt.ylabel('N10 Value', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Reverse the x-axis
plt.gca().invert_xaxis()


# Add legend, grid, and layout adjustments
plt.legend(loc='upper left', fontsize=14)
plt.title("Clustering Benchmark Results on Metabolomics Dataset", fontsize=16)
plt.grid(True)
plt.tight_layout()
plt.show()