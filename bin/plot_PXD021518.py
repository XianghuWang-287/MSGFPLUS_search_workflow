from matplotlib import pyplot as plt


falcon_x = [111,175,241,353,441]
falcon_y = [0.9772818529214622,0.9729616468781468,0.9677066403165359,0.9612795335878034,0.9529934196480674]

maracluster_x = [1,1,1,1,2]
maracluster_y = [0.987057036685737,0.9789106539318153,0.9675225271195903,0.9524003914105789,0.8947213600917293]

mscluster_x_MSRT = [35,58,115,176,393,644]
mscluster_y_MSRT = [0.9494786246887178,0.9378301487712926,0.9144581455999324,0.8873352002875002,0.8288441306715686,0.6028541649442926]

falcon_x_DB = [111,175,241,353,441]
falcon_y_DB = [0.985937922736661,0.9805808899527653,0.9742808085154929,0.9667904907315589,0.9564200625179712]


maracluster_x_DB = [1,1,1,1,2]
maracluster_y_DB = [0.9997941916304597, 0.9994511776812257,0.9987651497827579,0.9980333866910588,0.9939858220900983]

mscluster_x_DB = [35,58,115,176,393,644]
mscluster_y_DB = [0.9795228831780485,0.9722924587789452,0.9587011517115986,0.9406695874689827, 0.898622337045917,0.7904405341599882,]



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
plt.ylabel('N10 Value', fontsize=12)

# Reverse the x-axis
plt.gca().invert_xaxis()


# Add legend, grid, and layout adjustments
plt.legend(loc='lower right', fontsize=12)
plt.title("Clustering Benchmark Results on Joint Proteomics Dataset")
plt.grid(True)
plt.tight_layout()
plt.show()