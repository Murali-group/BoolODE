"""
 This script performs k-means and DBSCAN clustering
 algorithms on a provided ExpressionData.CSV file to estimate
 the number of steady states for a particular boolean model.

 @Author: Neel Patel (neelspatel999@vt.edu)
 @Author: Madeline Shaklee (mshaklee@umassd.edu)
 @Author: Amelia Whitehead (ameliafw@vt.edu)
"""
import numpy as np
import argparse
import pandas

parser = argparse.ArgumentParser(description='Takes the ExpressionData.csv filepath '
                                             'and various options to cluster and visualize the data.')

parser.add_argument('-f', '--ExpressionDataFilePath', required=True, type=str,
                    help='Specify the path of the ExpressionData file.')
parser.add_argument('-e', '--Elbow', action='store_true',
                    help='Input -e if you want to perform a k-means elbow visualization.')
parser.add_argument('-s', '--Silhouette', action='store_true',
                    help='Input -s if you want to perform a k-means silhouette visualization.')
parser.add_argument('-d', '--DBSCAN', action='store_true',
                    help='Input -d if you want to perform DBSCAN clustering to obtain a DBSCAN_ClusterIDs.csv file.')
parser.add_argument('-u', '--SilhouetteUpperBound', type=int,
                    help='Input -u if you want to specify the upper bound (inclusive) for the number of clusters to consider '
                         'in silhouette plots. The default is 11. You must use this option along with -s, '
                         'and you must provide a number greater than 2. NOTE: We suggest starting'
                         ' with low numbers (less than 12) because the image '
                         'will increase in size as you increase the number of silhouettes to generate')


elbow_value = parser.parse_args().Elbow
silhouette_value = parser.parse_args().Silhouette
dbscan_value = parser.parse_args().DBSCAN
silhouette_upperbound = parser.parse_args().SilhouetteUpperBound
expression_data_location = parser.parse_args().ExpressionDataFilePath

# Defining the output path to be the same path as
# ExpressionData for any output file generated
output_file_path = expression_data_location[0:expression_data_location.rindex('/') + 1]

""""
 Iteratively performs k-means clustering, testing k values 2 ≤ k ≤ upper_bound,
 and generates a silhouette plot of each.

 @Param:        upper_bound is the number of maximum number of clusters to
                include in the silhouette visualization. This can be specified
                by the user by -u or --SilhouetteUpperBound. The default
                is 11. upper_bound must be k ≥ 3
                
 @Precondition: The number of cells in the ExpressionData file 
                must be greater than upper_bound
                
 @Note:         The package only allows for k ≥ 2 analyses

 @Citation:    https://www.scikit-yb.org/en/latest/api/cluster/silhouette.html 
"""
def k_means_silhouette(upper_bound):
    from sklearn.cluster import KMeans
    from yellowbrick.cluster import SilhouetteVisualizer
    import matplotlib.pyplot as plt

    # ExpressionData is a matrix, where rows are genes and columns are cells
    expression_df = pandas.read_csv(expression_data_location, sep=',', index_col=0)

    # Performs transpose for later clustering on the samples.
    # Later, we cluster the cell columns of ExpressionData
    expression_data_transpose = np.array(expression_df.T.to_numpy().tolist())

    fig, axs = plt.subplots(upper_bound - 1, figsize=(15, upper_bound*5))
    plot_count = 0
    for cluster_number in range(2, upper_bound + 1):
        k_means = KMeans(n_clusters=cluster_number)

        visualizer = SilhouetteVisualizer(k_means, axs[plot_count])
        visualizer.fit(expression_data_transpose)
        axs[plot_count].set_title("Silhouette of " + str(cluster_number) + " Clusters")
        axs[plot_count].set_xlabel('silhouette coefficient values')
        axs[plot_count].set_ylabel('cluster label')
        plot_count = plot_count + 1

    visualizer.show(outpath=output_file_path + "silhouette_visualization.png")
    print("Silhouette analyses generated a silhouette_visualization.png file.")

""""
 Iteratively performs k-means clustering, testing k values 2 ≤ k ≤ 11,
 and generates an elbow plot.

 @Precondition: Number of cells in the ExpressionData file 
                must be greater than 11

 @Note:         The package only allows for k ≥ 2 analyses
                
 @Citation:     https://www.scikit-yb.org/en/latest/api/cluster/elbow.html
"""
def k_means_elbow():
    from sklearn.cluster import KMeans
    from yellowbrick.cluster import KElbowVisualizer

    # ExpressionData is a matrix, where rows are genes and columns are cells
    expression_df = pandas.read_csv(expression_data_location, sep=',', index_col=0)

    # Performs transpose for later clustering on the samples.
    # We cluster the cell columns of ExpressionData
    expression_data_transpose = np.array(expression_df.T.to_numpy().tolist())

    model = KMeans()

    # KElbowVisualizer has 3 different metrics
    # (distortion, silhouette, and calinski_harabasz)
    # The calinski_harabasz score computes the ratio of
    # dispersion between and within clusters. Distortion could also
    # be used for this method.
    visualizer = KElbowVisualizer(model, k=(2, 12), metric='calinski_harabasz', timings=False)

    visualizer.fit(expression_data_transpose)
    visualizer.show(outpath=output_file_path + "elbow_visualization.png")
    print("Elbow analysis generated an elbow_visualization.png file.")

"""
 Performs DBSCAN on the ExpressionData file and generates a DBSCAN_ClusterIDs.csv file, which specifies 
 which cell belongs to what cluster, according to DBSCAN. All noise points are grouped together into 
 a separate cluster for visualization purposes.
 
 @Precondition: The number of samples (cells) needs to be greater than two 
                times the number of genes
                
 @Note:         DBSCAN does not work well with datasets of varying density

 @Citations:    https://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html#sphx-glr-auto-examples-cluster-plot-dbscan-py
                https://medium.com/@tarammullin/dbscan-parameter-estimation-ff8330e3a3bd
"""
def dbscan_clustering():
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    from kneed import KneeLocator

    # ExpressionData is a matrix, where rows are genes and columns are cells
    expression_df = pandas.read_csv(expression_data_location, sep=',', index_col=0)

    # Performs transpose for later clustering on the samples.
    # We cluster the cell columns of ExpressionData
    expression_data_transpose = np.array(expression_df.T.to_numpy().tolist())


    # DBSCAN takes two parameters, MinSamples and Epsilon
    # MinSamples is described now, and Epsilon will be described later

    # The MinSamples parameter is the fewest number of samples (cells)
    # DBSCAN will put in a cluster. We assign MinSamples
    # to two times the number of genes in ExpressionData
    min_samples_num = len(expression_df.index)*2


    ############# Below is from Medium.com reference #############

    # Compute nearest neighbors and distances
    neighbors = NearestNeighbors(n_neighbors=min_samples_num)

    neighbors_fit = neighbors.fit(expression_data_transpose)
    distances, indices = neighbors_fit.kneighbors(expression_data_transpose)

    # Sort distances in ascending order
    distances = np.sort(distances, axis=0)

    ############# End of Medium.com reference #############


    # Find the average k-distances then plot to find the knee,
    # which is the point of maximum curvature
    y_values = compute_average_distances(distances)
    x_values = list(range(1,len(y_values) + 1))
    knee_locator = KneeLocator(x_values, y_values, S=1.0, curve='convex', direction='increasing')
    maximum_curvature_position = round(knee_locator.knee, 20)

    # The Epsilon parameter of DBSCAN represents the upper bound (exclusive)
    # distance between two points to be clustered together.
    epsilon = y_values[maximum_curvature_position - 1]


    ############# Below is from the scikit-learn reference #############

    # Performs DBSCAN with the calculated parameters
    db = DBSCAN(eps=epsilon, min_samples=min_samples_num).fit(expression_data_transpose)

    clusters_identifiers = db.fit_predict(expression_data_transpose)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print("DBSCAN result:")
    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d \n' % n_noise_)

    ############# End of scikit-learn reference #############

    # Adding the cluster labels to the dataframe
    expression_df.loc['cluster_labels', :] = clusters_identifiers

    # Write a new ClusterID file according to DBSCAN's clustering
    cluster_values = expression_df.loc['cluster_labels'].tolist()
    cluster_values = [int(number) for number in cluster_values]
    cell_names = expression_df.columns.tolist()
    dictionary = {'':cell_names, 'cl':cluster_values}
    cluster_df = pandas.DataFrame(dictionary)
    cluster_df.to_csv(output_file_path + 'DBSCAN_ClusterIDs.csv', index=False)

    print("DBSCAN analysis generated a DBSCAN_ClusterIDs.csv file.")

"""
 Takes a list of distance lists computed via k-nearest neighbors, and
 outputs a list where each element is the average distance of the distance 
 list occupied in the respective ordinal location of the input list
 
 @Param:        A list of lists, where each element contains a 
                list of distances computed via k-nearest neighbor.
 
 @Precondition: The distance_array has length greater than 0
 
 @Returns:      A list of distances, where each element is the 
                average distance of the distance list occupied 
                in the respective ordinal location of the input list
"""
def compute_average_distances(distance_array):
    denominator_of_average = len(distance_array[0])
    # Will be of length number of samples (cells) in the
    # ExpressionData, once full
    average_distances = []
    for distances in distance_array:
        distance_sum = 0
        for distance_value in distances:
            distance_sum = distance_sum + distance_value
        average_distance = distance_sum / denominator_of_average
        average_distances.append(average_distance)
    return average_distances

if __name__ == '__main__':
    if elbow_value:
        k_means_elbow()
    if silhouette_value:
        if silhouette_upperbound:
            k_means_silhouette(silhouette_upperbound)
        else:
            # The default value is 11
            k_means_silhouette(11)
    if dbscan_value:
        dbscan_clustering()
