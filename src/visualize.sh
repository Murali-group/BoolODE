echo "Starting time course visualization"
python ./timecourse-vis.py -i $1

echo "Starting dimensionality reduction visualization"
python ./dim-red-vis.py  -i $1

echo "Starting gene cluster visualization"
python ./cluster-vis.py -i $1

mv *.png ../

