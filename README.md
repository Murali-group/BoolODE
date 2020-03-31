# BoolODE
Git Repo for converting Boolean models to ODE models and performing stochastic simulations.

Find the documentation for BoolODE at [https://murali-group.github.io/Beeline/BoolODE.html](https://murali-group.github.io/Beeline/BoolODE.html).

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
