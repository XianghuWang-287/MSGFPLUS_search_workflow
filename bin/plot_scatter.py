from matplotlib import pyplot as plt

# Data
x_mscluster = [1743]
y_mscluster = [0.9495883807121245]
x_falcon = [2877]
y_falcon = [0.9642936116070135]
x_maracluster = [33]
y_maracluster = [0.9987753429790129]

# Create figure and plot
plt.figure(figsize=(10, 8))
plt.scatter(x_mscluster, y_mscluster, label='MSCluster', s=100)
plt.scatter(x_falcon, y_falcon, label='Falcon', s=100)
plt.scatter(x_maracluster, y_maracluster, label='MaRaCluster', s=100)

# Annotate each point with adjusted positions
annotations = [
    (x_mscluster, y_mscluster, 'MSCluster', (0, -25)),  # slightly further below
    (x_falcon, y_falcon, 'Falcon', (-20, -20)),  # shifted to the left
    (x_maracluster, y_maracluster, 'MaRaCluster', (34, -25))  # shifted to the right
]

for x, y, label, offset in annotations:
      # plot the points again to ensure they are above the annotations
    for i in range(len(x)):
        plt.annotate(f'({x[i]}, {y[i]:.3f})',
                     (x[i], y[i]),
                     textcoords="offset points",
                     xytext=offset,
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
