import numpy as np
import matplotlib.pyplot as plt
import hdbscan
from sklearn.datasets import make_moons
import matplotlib.cm as cm

# Create a dataset with two clusters of varying density
X, _ = make_moons(n_samples=300, noise=0.05, random_state=42)

# Initialize and fit HDBSCAN
clusterer = hdbscan.HDBSCAN(min_cluster_size=5, min_samples=3)
cluster_labels = clusterer.fit_predict(X)

# Print unique cluster labels for debugging
print("Unique cluster labels:", np.unique(cluster_labels))

# Generate a list of colors from the 'viridis' colormap
num_clusters = len(np.unique(cluster_labels))
colors = cm.get_cmap('viridis', num_clusters)

# Plot the data with cluster labels
plt.figure(figsize=(12, 6))

# Plot clusters
plt.subplot(1, 2, 1)
plt.title('HDBSCAN Clustering')
scatter = plt.scatter(X[:, 0], X[:, 1], c=cluster_labels, cmap='viridis', s=50)
plt.colorbar(scatter, label='Cluster Label')
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')

# Plot the condensed tree graph
plt.subplot(1, 2, 2)
plt.title('HDBSCAN Condensed Tree')
clusterer.condensed_tree_.plot(select_clusters=True, selection_palette=colors.colors)
plt.xlabel('Cluster Leaf Nodes')
plt.ylabel('Persistency')

# Save the plot as a file
plt.tight_layout()
plt.savefig('hdbscan_clusters.png')
plt.close()

print("Plots saved as 'hdbscan_clusters.png'")