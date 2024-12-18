import matplotlib.pyplot as plt

# Data
N10 = [53, 129, 180, 222, 277]
Completeness = [0.8438354712833923, 0.8782508581969166, 0.8977797710226826, 0.911072700517591, 0.9200205109252527]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(Completeness, N10, marker='o', linestyle='-', color='b')
plt.title('Correlation Analysis of Completeness and N10')
plt.xlabel('Completeness')
plt.ylabel('N10')
plt.grid(True)
plt.show()
