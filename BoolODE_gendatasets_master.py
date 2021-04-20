#BoolODE_gendatasets_master.py

import argparse
import sys
import os

pathtoBoolODE = '/home/cbuck016/BoolODE-0.1/src/'
outputdir = '/home/cbuck016/BoolODE-0.1/Test1/'

#def create_arg_parser():
#    parser = argparse.ArgumentParser(description='input pathway and output directory')
#    parser.add_argument('pathtoBoolODE', help='Path to BoolODE')
#    parser.add_argument('--output-dir', help='Path to output files')
#    return parser

#if __name__ == "__main__":
#    arg_parser = create_arg_parser()
#    parsed_args = arg_parser.parse_args(sys.argv[1:])
#    if os.path.exists(parsed_args.pathtoBoolODE):
#      print("File exists")

num_simulations = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
num_cells = [250, 500, 1000, 2000, 5000]
num_timesteps = [100, 200, 400, 800, 1600]
model_name = 'mCAD'

for cells in num_cells:
    for ts in num_timesteps:
        for i in num_simulations:
            experiment_name = model_name + '-ts-' + str(num_timesteps) + '-cells-' + str(num_cells) + '-sim-' + str(num_simulations)
            input_file_prefix = pathtoBoolODE + '/data/' + model_name 
            command = 'python3 BoolODE.py --path ' + input_file_prefix + '.txt' \
            + '--ics ' + input_file_prefix + '_ics.txt' \
            + '--max-time ' + str(num_timesteps) + ' --num-cells ' + str(num_cells) \
            + ' --do-parallel ' \
            + ' --outPrefix ' + outputdir + experiment_name + "/sim-" + str(num_simulations) + ".out" \
            + ' --sample-cells'
            print(command)
            os.system(command)
