echo "Starting time course visualization"
python src/timecourse-vis.py -i $1

echo "Starting dimensionality reduction visualization"
python src/dim-red-vis.py  -i $1

echo "Starting gene cluster visualization"
python src/cluster-vis.py -i $1


