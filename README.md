*NOTE: The BoolODE code is being refactored, and is currently being tested. Please use the V0.1 release for the version of the code used in the BEELINE publication*


# BoolODE
Git Repo for converting Boolean models to ODE models and performing stochastic simulations.

Find the documentation for BoolODE at [https://murali-group.github.io/Beeline/BoolODE.html](https://murali-group.github.io/Beeline/BoolODE.html).

## Installation

1) Clone the repository to get a copy of all BoolODE files on your personal machine.
In the directory where you want the BoolODE files, use the following command:
`git clone https://github.com/Murali-group/BoolODE.git`

2) Install the required software packages with the following command:
`pip install -r requirements.txt`

3) If you want to perform the optional post_processing of comparing DBSCAN with PyBoolNet (compare_pyboolnet:True), follow
the instructions to install the requirements for PyBoolNet here:
`https://pyboolnet.readthedocs.io/`

## Usage
`python boolode.py --config path/to/config.yaml`

## Configuration 
BoolODE reads user defined configurations from a YAML file. A sample config file is provided
in (config-files/example-config.yaml)[github.com//Murali-group/BoolODE/blob/master/config-files/example-config.yaml].

A minimal configuration would look like this:

```yaml
global_settings:
  model_dir: "data"
  output_dir: "Synthetic-H"
  do_simulations: True
  do_post_processing: False
  modeltype: 'hill'

jobs:
  - name: "name-of-run"
    model_definition: "boolean-model-definition"
    model_initial_conditions: "initial-condition-file"
    simulation_time: 9
    num_cells: 300
    do_parallel: True
    sample_cells: False
    
post_processing:
  - GenSamples:
    - sample_size: 200
      nDatasets: 10
      
  - DimRed:
    - perplexity: 200
    
  - Dropouts:
    - dropout: True
      sample_size: 100
      drop_cutoff: 0.5
      drop_prob: 0.5
  ```

## Inputs
BoolODE requires a tab separated file containing a Boolean model. The following is a sample input file
```
Gene	Rule
a	not(c)
c	not(b)
b	not(a)
```
The first column specifies the Gene, while the second column specifies the Boolean rule governing 
the expression of the gene. The rules should be formatted as follows: 
`(a_1 or a_2 ... ) and not (r_1 or r_2 ... )` where `a_`s are the activators, and `r_`s are the
repressors of a given gene.

## Outputs
BoolODE carries out as many SDE simulations as the number of cells requested. The trajectories of these simulations are stored in the `/simulations/` folder where they can be resampled. The simulation output relevant for use by GRN inference algorithms are the following:
1. `refNetwork.csv` - An edgelist with signs of interactions inferred from the model file.
2. `PseudoTime.csv` - A ground truth pseudotime file. BoolODE uses simulation time as a proxy for pseudotime. 
3. `ExpressionData.csv` - The table of gene expression values per 'cell'. For explanation of the format, see below.
4. `ClusterIds.csv` - A table assigning a cluster ID to each simulated trajectory by carrying out k-means clustering.

Additionally, BoolODE creates a `model.py` and a `parameters.txt` containing the the ODE model to be simulated, and the kinetic parameters values used to parameterize the ODE model.

The ExpressionData.csv file has rows corresponding to the genes, and columns corresponding to the timepoints in each experiment.  For example, `[E0_0,E0_10,E0_20,E1_0,E1_10,E1_20]` shows two experiments with 3 timepoints, at times 0,10,20 respectively.

## Overview of method
BoolODE is currently designed for transcription factor regulatory networks, though a protein interaction network can be specified by using the `species_type` option in the config file pointing to a tab separated file indicating the type of each variable, either `protein` or `gene`.

For each gene in the TF regulatory network, BoolODE creates two equations, one govering the regulation  of the gene (x), and one for its corresponding protein (p).
 
dx/dt = m f(X) - l\_x x
 
dp/dt = r x  - l\_p p

More more details, please see the supplementary material and Online Methods of the accompanying publication.

If you use BoolODE in your research, please cite:

Aditya Pratapa, Amogh Jalihal, Jeffrey Law, Aditya Bharadwaj, and T M Murali. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data, Nature Methods, (2020). https://doi.org/10.1038/s41592-019-0690-6

# Optional clustering techniques

## Perform_Clustering 

Performs clustering on an ExpressionData file with various cluster and visualization options. 
Clustering is used to predict the number of steady states in the Boolean Model expressed in the ExpressionData file. A brief explanation of the DBSCAN clustering technique is provided in the Overview section

## Usage

## Perform DBSCAN and create a DBSCAN_ClusterIDs.csv file

`python perform_clustering.py -f path/to/ExpressionData.csv -d`

## Perform k-means and creates an elbow_visualization.png file for clusters 2 to 11

`python perform_clustering.py -f path/to/ExpressionData.csv -e`

## Perform k-means silhouette and create a silhouette_visualization.png file for clusters 2 to 5
(Note: The default, without -u (upper bound) specified, is clusters 2 to 11)

`python perform_clustering.py -f path/to/ExpressionData.csv -s -u 5`

Note: any combination of -s, -e, and -d can be specified, for example:

`python perform_clustering.py -f path/to/ExpressionData.csv -d -e -s`


## Inputs

-f: The ExpressionData file path
Required
The ExpressionData file is expected to be the same format as generated by BoolODE

-d: DBSCAN clustering is requested
Not required

-e: k-means elbow is requested
Not required

-s: k-means silhouette is requested
Not required

-u: upper bound for silhouette plots, which must be greater than 2
Not required


## Outputs

-d: DBSCAN clustering is requested
The DBSCAN clustering outputs a DBSCAN_ClusterIDs.csv file, where each cell has its own cluster assigned by DBSCAN. Note: noise found by DBSCAN are put in a cluster for visualization purposes.

-e: k-means elbow is requested
The k-means elbow generates an elbow_visualization.png file with the elbow identified.

-s: k-means silhouette is requested
The k-means silhouette generates a silhouette_visualization.png file where the number of clusters is 2 until the number specified by -u. If no upper bound is specified by -u, the default upper bound is 11.

Note: All outputs are generated in the same directory as the provided ExpressionData.csv


## OUTPUT NOTES

## DBSCAN 

There is more work that needs to be done with DBSCAN, such as testing with lower simulation times and differing cell numbers. Additionally, DBSCAN is not suitable for datasets with varying density, and that is an analysis metric we want to include in the future (analyzing the density of the dataset and determining when DBSCAN is successful and what density works with particular models).

## ELBOW

The elbow method has one additional metric that can be implemented if the line defining the model type is changed from metric='calinski_harabasz' to metric='distortion’ in k_means_elbow(). We have set it to calinski_harabasz because that has yielded closer results to the expected value on average for the particular models that we tested.

## SILHOUETTE

The silhouette_visualization.png file needs to be visually analyzed to determine the appropriate number of clusters. The correct estimated cluster number, according to Silhouette, corresponds to the plot where there that has the least number of negative coefficient values and has the most uniformity in the thickness of the clusters.


## Overview

The perform_clustering script is designed to analyze ExpressionData files generated by BoolODE to statistically analyze the data to try and cluster it and predict the number of steady states that exist in the boolean model. It currently supports k-means elbow, k-means silhouette, and DBSCAN.

## DBSCAN Description

DBSCAN stands for Density-Based Spatial Clustering of Applications with Noise. It is useful for handling the data in ExpressionData files because of how it handles noise. DBSCAN works by looking at a point and a specified number of its neighbors at a certain distance to determine if they are similar enough to be grouped together in a cluster. In the end, DBSCAN outputs the estimated number of clusters that it found from the data. DBSCAN is used to estimate the number of steady states from simulated gene expression data from an ExpressionData file.
