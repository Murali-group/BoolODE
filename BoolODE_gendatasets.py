#generating datasets for Boolean models

path_to_boolode = "/home/cbuck016/BoolODE-0.1/"

maxtime = "8"
numcells = "1000"
numdatasets = "10"

output_dir="CuratedData/"


# Boolean model 1: mCAD
model_name = "mCAD"
output_name = "mCAD"
list_of_mCAD_sims = ["mCAD-sim-01-ts-800-cells-1000", "mCAD-sim-02-ts-800-cells-1000", "mCAD-sim-03-ts-800-cells-1000", "mCAD-sim-04-ts-800-cells-1000", "mCAD-sim-05-ts-800-cells-1000", "mCAD-sim-06-ts-800-cells-1000", "mCAD-sim-07-ts-800-cells-1000", "mCAD-sim-08-ts-800-cells-1000", "mCAD-sim-09-ts-800-cells-1000", "mCAD-sim-10-ts-800-cells-1000"]
for mCAD_sim in list_of_mCAD_sims:
	echo "Simulating "$model_name
	python3 $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name".txt" \
						--ics $path_to_boolode/data/$model_name"_ics.txt" \
						--max-time $maxtime --num-cells $numcells \
						--do-parallel \
						--outPrefix $output_dir$output_name"/" \
						--sample-cells 
	

# Boolean model 2: VSC
model_name = "VSC"
output_name = "VSC"
list_of_VSC_sims = ["VSC-sim-01-ts-800-cells-1000", "VSC-sim-02-ts-800-cells-1000", "VSC-sim-03-ts-800-cells-1000", "VSC-sim-04-ts-800-cells-1000", "VSC-sim-05-ts-800-cells-1000", "VSC-sim-06-ts-800-cells-1000", "VSC-sim-07-ts-800-cells-1000", "VSC-sim-08-ts-800-cells-1000", "VSC-sim-09-ts-800-cells-1000", "VSC-sim-10-ts-800-cells-1000"]
for VSC_sim in list_of_VSC_sims:
	echo "Simulating "$model_name
	python3 $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name".txt" \
						--max-time $maxtime --num-cells $numcells \
						--do-parallel \
						--outPrefix $output_dir$output_name"/" \
						--sample-cells 


# Boolean model 3: HSC
model_name = "HSC"
output_name = "HSC"
list_of_HSC_sims = ["HSC-sim-01-ts-800-cells-1000", "HSC-sim-02-ts-800-cells-1000", "HSC-sim-03-ts-800-cells-1000", "HSC-sim-04-ts-800-cells-1000", "HSC-sim-05-ts-800-cells-1000", "HSC-sim-06-ts-800-cells-1000", "HSC-sim-07-ts-800-cells-1000", "HSC-sim-08-ts-800-cells-1000", "HSC-sim-09-ts-800-cells-1000", "HSC-sim-10-ts-800-cells-1000"]
for HSC_sim in list_of_HSC_sims:
	echo "Simulating "$model_name
	python3 $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name".txt" \
						--ics $path_to_boolode/data/$model_name"_ics.txt" \
						--max-time $maxtime --num-cells $numcells \
						--do-parallel \
						--outPrefix $output_dir$output_name"/" \
						--sample-cells 


# Boolean model 4: GSD
model_name = "GSD"
output_name = "GSD"
list_of_GSD_sims = ["GSD-sim-01-ts-800-cells-1000", "GSD-sim-02-ts-800-cells-1000", "GSD-sim-03-ts-800-cells-1000", "GSD-sim-04-ts-800-cells-1000", "GSD-sim-05-ts-800-cells-1000", "GSD-sim-06-ts-800-cells-1000", "GSD-sim-07-ts-800-cells-1000", "GSD-sim-08-ts-800-cells-1000", "GSD-sim-09-ts-800-cells-1000", "GSD-sim-10-ts-800-cells-1000"]
for GSD_sim in list_of_GSD_sims:
	echo "Simulating "$model_name
	python3 $path_to_boolode/src/BoolODE.py --path $path_to_boolode/data/$model_name".txt" \
						--ics $path_to_boolode/data/$model_name"_ics.txt" \
						--max-time $maxtime --num-cells $numcells \
						--do-parallel \
						--outPrefix $output_dir$output_name"/" \
						--sample-cells 



