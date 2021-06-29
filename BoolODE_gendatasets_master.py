#BoolODE_gendatasets_master.py

import argparse
import sys
import os

pathtoBoolODE = '/home/cbuck016/BoolODE-0.1/'
outputdir = '/home/cbuck016/BoolODE-0.1/GSD-sims/'

num_simulations = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
num_cells = [250, 500, 1000, 2000, 5000]
num_timesteps = [1, 2, 4, 8, 16]
model_name = 'GSD'

for cells in num_cells:
    for ts in num_timesteps:
        for i in num_simulations:
            experiment_name = model_name + '-ts-' + str(ts) + '00' + '-cells-' + str(cells) + '-sim-' + str(i)
            input_file_prefix = pathtoBoolODE + 'data/' + model_name 
            command = 'python3 ' + pathtoBoolODE + 'src/BoolODE.py --path ' + input_file_prefix + '.txt ' \
            + ' --ics ' + input_file_prefix + '_ics.txt ' \
            + ' --max-time ' + str(ts) + ' --num-cells ' + str(cells) \
            + ' --do-parallel ' \
            + ' --outPrefix ' + outputdir + experiment_name + "/sim-" + str(i) + "-" \
            + ' --sample-cells'
            print(command)
            os.system(command)
