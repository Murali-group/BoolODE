# Set the path to your clone of the BoolODE repository here
# The current options will perform stochastic simulations
# using the default kinetic parameters.
# If you would like to sample parameters, use the following options:
#  --sample-par --std 0.05 -i
# The above will specify a standard deviation of 5% of the default parameter
# value, and the -i option will set all parameters to a single sampled value.

path_to_boolode="../"
NUMCELLS="2000"

# Synthetic Networks
# Linear
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-linear.txt\
       --ics $path_to_boolode/data/dyn-linear_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --outPrefix "Synthetic/dyn-LI/"

# Linear Long
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-linear-long.txt\
       --ics $path_to_boolode/data/dyn-linear-long_ics.txt --max-time 15 --num-cells $NUMCELLS\
       --do-parallel\
       --outPrefix "Synthetic/dyn-LL/"
       
# Cycle
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-cycle.txt\
       --ics $path_to_boolode/data/dyn-cycle_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --outPrefix "Synthetic/dyn-CY/"
       
# Bifurcating
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-bifurcating.txt\
       --ics $path_to_boolode/data/dyn-bifurcating_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "Synthetic/dyn-BF/"

# Bifurcating-converging
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-bifurcating-converging.txt\
       --ics $path_to_boolode/data/dyn-bifurcating-converging_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "Synthetic/dyn-BFC/"

# Trifurcating
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/dyn-trifurcating.txt\
       --ics $path_to_boolode/data/dyn-trifurcating_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 3\
       --outPrefix "Synthetic/dyn-TF/"



# Curated models
# mCAD
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/mCAD.txt\
       --ics $path_to_boolode/data/mCAD_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "Curated/mCAD/"

# VSC
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/VSC.txt\
       --ics $path_to_boolode/data/VSC_ics.txt --max-time 5 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 5\
       --outPrefix "Curated/VSC/"

# HSC
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/HSC.txt\
       --ics $path_to_boolode/data/HSC_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 4\
       --outPrefix "Curated/HSC/"

# GSD
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/GSD.txt\
       --ics $path_to_boolode/data/GSD_ics.txt --max-time 8 --num-cells $NUMCELLS\
       --do-parallel\
       --nClusters 2\
       --outPrefix "Curated/GSD/"
       
