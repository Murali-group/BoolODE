"""
 This script writes config files based on
 a given range of cells and range of time values for a given model.
 It creates a directory for each cell and time combination.

 Invoke this file as follows:
 python generate_config_files.py path_of_input_file.txt
 The example_gen_config_input.txt file is the format expected.

 @Author: Neel Patel (neelspatel999@vt.edu)
 @Author: Amelia Whitehead (ameliafw@vt.edu)
"""

import sys
import os

base_directory_path = ""
model_name = ""
model_definition = ""
model_initial_conditions = ""
cells_list = []
time_list = []
CompareSteadyStates = False
perform_PyBoolNet = False

# NOTE: Right now, we are only considering CompareSteadyStates
# by itself or along with perform_PyBoolNet as
# post_processing requests. We do not look for any other
# post_processing requests at the moment.
do_post_processing = False
def main():
    """
     Main caller function of the generate_config_files.py
    """

    for input_file in sys.argv[1:]:
        file = open(input_file, "r")
        lines = file.readlines()

        # Get information from input file. #
        ###########################################
        for line in lines:
            line = line.strip()
            colon_index = line.find(":")
            if "base_directory_path:" in line:
                global base_directory_path
                base_directory_path += line[colon_index+1:].strip()
            elif "name:" in line:
                global model_name
                model_name += line[colon_index+1:].strip()
            elif "model_definition:" in line:
                global model_definition
                model_definition += line[colon_index+1:].strip()
            elif "model_initial_conditions:" in line:
                global model_initial_conditions
                model_initial_conditions += line[colon_index+1:].strip()
            elif "simulation_times:" in line:
                times = line[colon_index+1:]
                # Split the string up by commas
                global time_list
                time_list = times.split(", ")
            elif "num_cells:" in line:
                cells = line[colon_index+1:]
                # Split the string up by commas
                global cells_list
                cells_list = cells.split(", ")
            elif "CompareSteadyStates:" in line:
                if "True" in line:
                    global CompareSteadyStates
                    CompareSteadyStates = True
            elif "perform_PyBoolNet:" in line:
                if "True" in line:
                    global perform_PyBoolNet
                    perform_PyBoolNet = True
        ###########################################
        generate_requested_config_files()

def generate_requested_config_files():
    """
     Creates a directory named the model name provided.
     Then creates subdirectories for each cell and time combination.
     Writes a config file in each directory created.
    """

    # Create the directory where all the files
    # will be located
    base_path = os.path.join(base_directory_path, model_name.strip('"'))
    os.mkdir(base_path)
    for cell_number in cells_list:
        # Create a directory for that cell number
        cell_number = cell_number.strip()
        cell_name = "Cells="
        cell_name += cell_number
        cell_specific_path = os.path.join(base_path, cell_name)
        os.mkdir(cell_specific_path)
        for time_amount in time_list:
            time_amount = time_amount.strip()
            time_point = "Time="
            time_point += time_amount
            time_specific_path = os.path.join(cell_specific_path, time_point)
            os.mkdir(time_specific_path)

            ## Write all the information in a config file.
            config_file_path = os.path.join(time_specific_path, "config.yaml")
            config_file = open(config_file_path, "w")
            config_file.write("global_settings:\n")
            config_file.write("  model_dir: data\n")
            info_to_write1 = "  output_dir: "
            info_to_write1 += time_specific_path
            info_to_write1 += "\n"
            config_file.write(info_to_write1)
            config_file.write("  do_simulations: True\n")
            if CompareSteadyStates:
                config_file.write("  do_post_processing: True\n")
            else:
                config_file.write("  do_post_processing: False\n")
            config_file.write("  modeltype: 'hill'\n")
            config_file.write("jobs:\n")
            info_to_write = "    - name: "
            info_to_write += model_name
            info_to_write += "\n"
            config_file.write(info_to_write)
            info_to_write = "      model_definition: "
            info_to_write += model_definition
            info_to_write += "\n"
            config_file.write(info_to_write)
            if model_initial_conditions != "" || model_initial_conditions != "None":
                info_to_write = "      model_initial_conditions: "
                info_to_write += model_initial_conditions
                info_to_write += "\n"
                config_file.write(info_to_write)
            info_to_write = "      simulation_time: "
            info_to_write += time_amount
            info_to_write += "\n"
            config_file.write(info_to_write)
            info_to_write = "      num_cells: "
            info_to_write += cell_number
            info_to_write += "\n"
            config_file.write(info_to_write)
            config_file.write("      do_parallel: False\n")
            config_file.write("      sample_cells: False\n")
            config_file.write("post_processing:\n")
            if CompareSteadyStates:
                config_file.write("  CompareSteadyStates:\n")
                info_to_write = "    - outputPath: "
                info_to_write += time_specific_path
                info_to_write += "\n"
                config_file.write(info_to_write)
                if perform_PyBoolNet:
                    config_file.write("      perform_PyBoolNet: True\n")
                else:
                    config_file.write("      perform_PyBoolNet: False\n")

if __name__ == "__main__":
    """
     Calls the main method
    """
    main()
