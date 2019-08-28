# Set the path to your clone of the BoolODE repository here
# The current options will perform stochastic simulations
# using the default kinetic parameters.
# If you would like to sample parameters, use the following options:
#  --sample-par --std 0.05 -i
# The above will specify a standard deviation of 5% of the default parameter
# value, and the -i option will set all parameters to a single sampled value.

path_to_boolode=".."
numcells="2000"

# Synthetic Networks
output_dir="Synthetic/"

# Linear
output_name="dyn-LI"
echo "Simulating "$model_name
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 5 --num-cells $numcells\
       --do-parallel\
       --outPrefix $output_dir$output_name"/"

# Linear Long
output_name="dyn-LL"
model_name="dyn-linear-long"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 15 --num-cells $numcells\
       --do-parallel\
       --outPrefix $output_dir$output_name"/"
       
# Cycle
output_name="dyn-CY"
model_name="dyn-cycle"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 8 --num-cells $numcells\
       --do-parallel\
       --outPrefix $output_dir$output_name"/"
       
# Bifurcating
output_name="dyn-BF"
model_name="dyn-bifurcating"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 5 --num-cells $numcells\
       --do-parallel\
       --nClusters 2\
       --outPrefix $output_dir$output_name"/"

# Bifurcating-converging
model_name="dyn-bifurcating-converging"
output_name="dyn-BFC"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 8 --num-cells $numcells\
       --do-parallel\
       --nClusters 2\
       --outPrefix $output_dir$output_name"/"

# Trifurcating
model_name="dyn-trifurcating"
output_name="dyn-TF"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 8 --num-cells $numcells\
       --do-parallel\
       --nClusters 3\
       --outPrefix $output_dir$output_name"/"


# # Curated models
output_dir="Curated/"

# mCAD
output_name="mCAD"
model_name="mCAD"
numclusters="2"
echo "Simulating "$model_name
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 5 --num-cells $numcells\
       --do-parallel\
       --nClusters $numclusters\
       --outPrefix $output_dir$output_name"/"
       
echo "Running slingshot"
numclusters=3
python $path_to_boolode/scripts/runSlingshot.py --expr $output_dir$output_name"/ExpressionData.csv"\
       --pseudo $output_dir$output_name"/PseudoTime.csv"\
       --outPrefix $output_dir$output_name\
       --nClusters=$numclusters -r 200


echo "Generating dropouts for q=50 for "$model_name
python $path_to_boolode/scripts/genDropouts.py  -e $output_dir$output_name"/ExpressionData.csv"\
       -p $output_dir$output_name"/PseudoTime.csv"\
       -r $output_dir$output_name"/refNetwork.csv"\
       -n $numcells -d --drop-cutoff 0.5 --drop-prob 0.5 -i 1\
       --outPrefix $output_dir$output_name"/"$output_name"-50/"


# VSC
model_name="VSC"
output_name="VSC"
numclusters="5"
echo "Simulating "$model_name
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 5 --num-cells $numcells\
       --do-parallel\
       --nClusters $numclusters\
       --outPrefix $output_dir$output_name"/"

echo "Running slingshot"
numclusters=6
python $path_to_boolode/scripts/runSlingshot.py --expr $output_dir$output_name"/ExpressionData.csv"\
       --pseudo $output_dir$output_name"/PseudoTime.csv"\
       --outPrefix $output_dir$output_name\
       --nClusters=$numclusters -r 200
       
echo "Generating dropouts for q=50 for "$model_name
python $path_to_boolode/scripts/genDropouts.py  -e $output_dir$output_name"/ExpressionData.csv"\
       -p $output_dir$output_name"/PseudoTime.csv"\
       -r $output_dir$output_name"/refNetwork.csv"\
       -n $numcells -d --drop-cutoff 0.5 --drop-prob 0.5 -i 1\
       --outPrefix $output_dir$output_name"/"$output_name"-50/"

# HSC
model_name="HSC"
output_name="HSC"
numclusters="4"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/HSC.txt\
       --ics $path_to_boolode/data/$model_name"_ics.txt" --max-time 8 --num-cells $numcells\
       --do-parallel\
       --nClusters $numclusters\
       --outPrefix $output_dir$output_name"/"
       
echo "Running slingshot"
numclusters=5
python $path_to_boolode/scripts/runSlingshot.py --expr $output_dir$output_name"/ExpressionData.csv"\
       --pseudo $output_dir$output_name"/PseudoTime.csv"\
       --outPrefix $output_dir$output_name\
       --nClusters=$numclusters -r 200
       
echo "Generating dropouts for q=50 for "$model_name
python $path_to_boolode/scripts/genDropouts.py  -e $output_dir$output_name"/ExpressionData.csv"\
       -p $output_dir$output_name"/PseudoTime.csv"\
       -r $output_dir$output_name"/refNetwork.csv"\
       -n $numcells -d --drop-cutoff 0.5 --drop-prob 0.5 -i 1\
       --outPrefix $output_dir$output_name"/"$output_name"-50/"



# GSD
output_name="GSD"
model_name="GSD"
numclusters="2"
python $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name.txt\
       --max-time 8 --num-cells $numcells\
       --do-parallel\
       --nClusters $numclusters\
       --outPrefix $output_dir$output_name"/"

echo "Running slingshot"
numclusters=3
python $path_to_boolode/scripts/runSlingshot.py --expr $output_dir$output_name"/ExpressionData.csv"\
       --pseudo $output_dir$output_name"/PseudoTime.csv"\
       --outPrefix $output_dir$output_name\
       --nClusters=$numclusters -r 200
       
echo "Generating dropouts for q=50 for "$model_name
python $path_to_boolode/scripts/genDropouts.py  -e $output_dir$output_name"/ExpressionData.csv"\
       -p $output_dir$output_name"/PseudoTime.csv"\
       -r $output_dir$output_name"/refNetwork.csv"\
       -n $numcells -d --drop-cutoff 0.5 --drop-prob 0.5 -i 1\
       --outPrefix $output_dir$output_name"/"$output_name"-50/"

