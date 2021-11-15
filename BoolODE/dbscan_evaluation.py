"""
 This script DBSCAN clustering
 on a provided ExpressionData.CSV file to estimate
 the number of steady states for a particular boolean model.

 @Author: Neel Patel (neelspatel999@vt.edu)
 @Author: Madeline Shaklee (mshaklee@umassd.edu)
 @Author: Amelia Whitehead (ameliafw@vt.edu)
"""

"""
 Performs DBSCAN on the ExpressionData file and generates a DBSCAN_ClusterIDs.csv file, which specifies 
 which cell belongs to what cluster, according to DBSCAN. All noise points are grouped together into 
 a separate cluster for visualization purposes.

 @Precondition: The number of samples (cells) needs to be greater than two 
                times the number of genes

 @Note:         DBSCAN does not work well with datasets of varying density
 
 @Params:       expression_data_location: the full file path of ExpressionData
                output_file_path: the file path where output files should go
                compare_with_pyboolnet: 1 if compare output with PyBoolNet, 0 otherwise
                model_path: States the path of the model. Needed for PyBoolNet to analyze steady states

 @Citations:    https://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html#sphx-glr-auto-examples-cluster-plot-dbscan-py
                https://medium.com/@tarammullin/dbscan-parameter-estimation-ff8330e3a3bd  
"""


def dbscan_clustering(expression_data_location, output_file_path, compare_with_pyboolnet, model_path):
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    from kneed import KneeLocator
    import numpy as np
    import pandas

    # ExpressionData is a matrix, where rows are genes and columns are cells
    expression_df = pandas.read_csv(expression_data_location, sep=',', index_col=0)

    # Performs transpose for later clustering on the samples.
    # We cluster the cell columns of ExpressionData
    expression_data_transpose = np.array(expression_df.T.to_numpy().tolist())
    # ExpressionData is a matrix, where rows are genes and columns are cells
    # expression_df = pandas.read_csv(expression_data_location, sep=',', index_col=0)

    # Performs transpose for later clustering on the samples.
    # We cluster the cell columns of ExpressionData
    # expression_data_transpose = np.array(data_frame.tolist())

    # DBSCAN takes two parameters, MinSamples and Epsilon
    # MinSamples is described now, and Epsilon will be described later

    # The MinSamples parameter is the fewest number of samples (cells)
    # DBSCAN will put in a cluster. We assign MinSamples
    # to two times the number of genes in ExpressionData
    min_samples_num = len(expression_df.index) * 2

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
    # y_values = compute_average_distances(distances)
    #############################################################
    denominator_of_average = len(distances[0])
    # Will be of length number of samples (cells) in the
    # ExpressionData, once full

    average_distances = []
    for distances_element in distances:
        distance_sum = 0
        for distance_value in distances_element:
            distance_sum = distance_sum + distance_value
        average_distance = distance_sum / denominator_of_average
        average_distances.append(average_distance)
    y_values = average_distances
    ############################################
    x_values = list(range(1, len(y_values) + 1))
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
    dictionary = {'': cell_names, 'cl': cluster_values}
    cluster_df = pandas.DataFrame(dictionary)
    cluster_df.to_csv(output_file_path + 'DBSCAN_ClusterIDs.csv', index=False)

    print("DBSCAN analysis generated a DBSCAN_ClusterIDs.csv file.")
    pyboolnet_result = -1
    if compare_with_pyboolnet:
        # Note: not the correct way to express the file path.
        pyboolnet_result = pyBoolNet_comparison(0, False, model_path)

    # Only do these print statements if PyBoolNet was ran.
    if pyboolnet_result != -1:
        if n_clusters_ == pyboolnet_result:
            print("\nDBSCAN and PyBoolNet results match and found " + str(n_clusters_) + "steady states\n")
        else:
            print("\nDBSCAN and PyBoolNet outputs did not match. \nDBSCAN: " + str(n_clusters_) + ", PyBoolNet: " + str(
                pyboolnet_result) + "\n")


# Runs PyBoolNet
def pyBoolNet_comparison(bnet_file, ics_present, model_path):
    #### Still has initial conditions section to be added
    print("Starting PyBoolNet...")

    import PyBoolNet
    import pandas as pd
    import os

    input_file = open(model_path, "r")
    temp_bnet_file = open(model_path + "_temporary_bnet.txt", "x")
    for line in input_file:
        # If it is not the first line
        if "".join(line.split()).lower() != "generule":
            # Convert the format and remove tabs.
            py_formatted_line = line.replace('\t', ',', 1)
            py_formatted_line = py_formatted_line.replace("and", "&")
            py_formatted_line = py_formatted_line.replace("or", "|")
            py_formatted_line = py_formatted_line.replace("not", "!")
            py_formatted_line = py_formatted_line.replace('\t', ' ')
            temp_bnet_file.write(py_formatted_line)
            temp_bnet_file.write('\n')

    input_file.close()
    temp_bnet_file.close()

    primes = PyBoolNet.FileExchange.bnet2primes(model_path + "_temporary_bnet.txt")
    steady_states = PyBoolNet.AspSolver.steady_states(primes)
    df = pd.DataFrame(steady_states)

    print("PyBoolNet found " + str(len(df)) + " steady state/s total")

    # Delete the temporary file.
    os.remove(model_path + "_temporary_bnet.txt")

    # Initial conditions and reachability from a given initial condition is TBD
    # We currently always set this to always false. Still need to work on this section
    unreachable_steady_states = 0
    if ics_present:
        current_goal_state = ""
        ic = "1--00"
        list_of_steady_states = df.values.tolist()
        # Find the number of unreachable_steady_states
        for steady_state in list_of_steady_states:
            current_goal_state = ','.join(str(v) for v in steady_state)
            current_goal_state = current_goal_state.replace(',', '')
            current_goal_state = current_goal_state.strip()
            path = PyBoolNet.StateTransitionGraphs.best_first_reachability(primes, InitialSpace=ic,
                                                                           GoalSpace=current_goal_state)
            if path:
                # for x in path:
                #     print(x)
                print("A path was found for steady state: " + current_goal_state)
            else:
                print("No path was found for steady state: " + current_goal_state)

    return len(df) - unreachable_steady_states
