#Silhouette Model example code below retrieved from scitkit-learn, currently editing to fit/run on BoolODE expression data

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import matplotlib.style as style

#Read ExpressionData.csv file
data = pd.read_csv (r'--FilePathtoExpressionData.csv--')
X = data.drop(data.columns[0], axis=1).transpose()
print(X)

range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10]
silhouette_avg_n_clusters = []

for n_clusters in range_n_clusters:
    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    clusterer = KMeans(n_clusters=n_clusters, random_state=42)
    cluster_labels = clusterer.fit_predict(X)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    silhouette_avg = silhouette_score(X, cluster_labels)
    print("For n_clusters =", n_clusters,
          "The average silhouette_score is :", silhouette_avg)

    silhouette_avg_n_clusters.append(silhouette_avg)
    
best_avg_silhouette_value = max(silhouette_avg_n_clusters)
best_num_cluster = silhouette_avg_n_clusters.index(best_avg_silhouette_value)

print("The best average silhouette score is: ", best_avg_silhouette_value)
