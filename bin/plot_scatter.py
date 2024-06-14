from matplotlib import pyplot as plt

# Data
# x_mscluster = [1743]
# y_mscluster = [0.9495883807121245]
# x_falcon = [2877]
# y_falcon = [0.9642936116070135]
# x_maracluster = [33]
# y_maracluster = [0.9987753429790129]
x_mscluster = [586]
y_mscluster = [0.8404377977378951]
x_falcon = [96]
y_falcon = [0.9780419296350973]
x_maracluster = [82]
y_maracluster = [0.9782067017271998]

y_mscluster_DB = [0.8446518259382786]
y_falcon_DB = [0.9842243711050124]
y_maracluster_DB = [0.9837496485947217]


# Create figure and plot
plt.figure(figsize=(10, 8))
plt.scatter(x_mscluster, y_mscluster, label='MSCluster_MSRT', s=100)
plt.scatter(x_falcon, y_falcon, label='Falcon_MSRT', s=100)
plt.scatter(x_maracluster, y_maracluster, label='MaRaCluster_MSRT', s=100)
plt.scatter(x_mscluster, y_mscluster_DB, label='MSCluster_DB', s=100)
plt.scatter(x_falcon, y_falcon_DB, label='Falcon_DB', s=100)
plt.scatter(x_maracluster, y_maracluster_DB, label='MaRaCluster_DB', s=100)

# Annotate each point with adjusted positions
annotations = [
    (x_mscluster, y_mscluster, 'MSCluster_MSRT'),  # slightly further below
    (x_falcon, y_falcon, 'Falcon_MSRT'),  # shifted to the left
    (x_maracluster, y_maracluster, 'MaRaCluster_MSRT'),  # shifted to the right
    (x_mscluster, y_mscluster_DB, 'MSCluster_DB'),  # slightly further below
    (x_falcon, y_falcon_DB, 'Falcon_DB'),  # shifted to the left
    (x_maracluster, y_maracluster_DB, 'MaRaCluster_DB')  # shifted to the right
]

for x, y, label in annotations:
      # plot the points again to ensure they are above the annotations
    for i in range(len(x)):
        plt.annotate(f'({x[i]}, {y[i]:.3f})',
                     (x[i], y[i]),
                     textcoords="offset points",
                     xytext=(80,-10),
                     ha='center',
                     fontsize=12,
                     bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))

# Adjust plot boundaries and set the x-axis to start from 0
plt.gca().margins(0.2)
plt.xlim(left=0)

# Labels and grid
plt.xlabel('N50 Value', fontsize=20)
plt.ylabel('Weighted Average Purity', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=15)
plt.grid(True)

# Show plot
plt.show()
