# BoolODE
Git Repo for simulating Boolean Models

# Usage
`python src/BoolODE.py --path=data/variables.txt --max-time=200 --num-timepoints=10 --num-experiments=3`

# Output
Produces two files by default, `stoch_experiment.txt` and `ode_experient.txt`. Each output
file has rows corresponding to the genes, and columns corresponding to the timepoints in each experiment.
For example, `[E0_0,E0_10,E0_20,E1_0,E1_10,E1_20]` shows two experiments with
3 timepoints, at times 0,10,20 respectively.
