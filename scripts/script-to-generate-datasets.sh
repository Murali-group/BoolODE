# Set the path to your clone of the BoolODE repository here
# The current options will perform stochastic simulations
# using the default kinetic parameters.
# If you would like to sample parameters, use the following options:
#  --sample-par --std 0.05 -i
# The above will specify a standard deviation of 5% of the default parameter
# value, and the -i option will set all parameters to a single sampled value.

path_to_boolode="../"
NUMCELLS="5000"

# Dyn- models
# linear
## simulate and store trajectories
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-linear.txt\
       --ics $path_to_boolode/data/dyn-linear_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --outPrefix "generated-datasets/dyn-linear/"
## Carry out dropouts
python $path_to_boolode/src/genDropouts.py  -e "benchmarking/Simulated-samples/"$MODELNAME"/"$MODELNAME"-"$NUMCELLS"-"$simnum"/ExpressionData.csv"\
       -p "benchmarking/Simulated-samples/"$MODELNAME"/"$MODELNAME"-"$NUMCELLS"-"$simnum"/PseudoTime.csv"\
       -r "benchmarking/Simulated-samples/"$MODELNAME"/"$MODELNAME"-"$NUMCELLS"-"$simnum"/refNetwork.csv"\
       -n 2000 -d --drop-cutoff 0.7 --drop-prob 0.7 -i $simnum\
       --outPrefix "benchmarking/Simulated-dropouts/"$MODELNAME"/"$MODELNAME
# linear long
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-linear-long.txt\
       --ics $path_to_boolode/data/dyn-linear-long_ics.txt --max-time 15 --num-cells $NUMCELLS\
       --do-parallel\
       --outPrefix "generated-datasets/dyn-linear-long/"

# bifurcating
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-bifurcating.txt\
       --ics $path_to_boolode/data/dyn-bifurcating_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "generated-datasets/dyn-bifurcating/"

# bifurcating converging
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-bifurcating-converging.txt\
       --ics $path_to_boolode/data/dyn-bifurcating-converging_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "generated-datasets/dyn-bifurcating-converging/"

# trifurcating
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-trifurcating.txt\
       --ics $path_to_boolode/data/dyn-trifurcating_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 3\
       --outPrefix "generated-datasets/dyn-trifurcating/"

# cycle
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-cycle.txt\
       --ics $path_to_boolode/data/dyn-cycle_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --outPrefix "generated-datasets/dyn-cycle/"

# Boolean models
# mCAD
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/mCAD.txt\
       --ics $path_to_boolode/data/mCAD_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "generated-datasets/mCAD/"

# VSC
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/VSC.txt\
       --ics $path_to_boolode/data/VSC_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 5\
       --outPrefix "generated-datasets/VSC/"

# HSC
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/HSC.txt\
       --ics $path_to_boolode/data/HSC_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 4\
       --outPrefix "generated-datasets/HSC/"

# GSD
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/GSD.txt\
       --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "generated-datasets/HSC/"
       
