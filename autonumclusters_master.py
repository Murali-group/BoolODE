#Silhouette Method for use with output results from BoolODE datasets -- IN PROGRESS

import os

pathtoBoolOutFile = '/home/cbuck016/BoolODE-0.1/mCAD-sims/'

num_simulations = [1]
num_cells = [250]
num_timesteps = [1]
#num_simulations = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#num_cells = [250, 500, 1000, 2000, 5000]
#num_timesteps = [1, 2, 4, 8, 16]
model_name = 'mCAD'

for cells in num_cells:
  for ts in num_timesteps:
    for i in num_simulations:
        experiment_name = model_name + '-ts-' + str(ts) + '00' + '-cells-' + str(cells) + '-sim-' + str(i)
        command = 'python3 ' + 'autonumclusters.py --expression-file ' + pathtoBoolOutFile + experiment_name \
        + '/' + 'sim-' + str(i) + '-ExpressionData.csv ' + '--outPrefix ' + pathtoBoolOutFile + experiment_name + "/" 
        print(command)
        os.system(command)
