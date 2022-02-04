__author__ = 'Jon Mallen'

# This is the common binarization method used by genSS.py and by extension, genVis.py. The expression data is sorted
# into an ascending list and clustered into 2 clusters using k-means. The average of the centroids of the two clusters
# is taken to be the binarization threshold for all of the data.

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans


def binarize_data(dataframe):
    gene_values = dataframe.values.tolist()
    all_values = [item for sublist in gene_values for item in sublist]
    all_values.sort()
    value_coordinate_list = [[a, all_values[a]] for a in range(len(all_values))]
    values_array = np.array(value_coordinate_list)
    expression_k_means = KMeans(n_clusters=2).fit(values_array)
    threshold = np.float32(round(sum(expression_k_means.cluster_centers_[:, 1]) / 2, 3))
    if threshold > 1:
        threshold = 1
    binarized_dataframe = pd.DataFrame()
    for c in dataframe.columns:
        binarized_dataframe[c] = (dataframe[c] >= threshold).astype(int)
    return binarized_dataframe
