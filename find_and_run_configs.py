import argparse
import BoolODE as bo
import os

def get_input_info():
    """
     Gets the root_directory path from the user.
    """

    parser = argparse.ArgumentParser(
        description='Run BoolODE to generate synthetic scRNAseq data for all'
                    'config files located within the given directory and '
                    'its subdirectories.\n '
                    'Specify the main directory where all config files are'
                    'contained. This is to be used in conjunction with '
                    'generate_config_files.py.')

    parser.add_argument('--root_directory', default="None",
                        help='Provide the root directory path'
                             'created by generate_config_files.py. '
                             'where its subdirectories contain config files.')

    return parser.parse_args().root_directory


def main():
    """
     Main caller function of the find_and_run_configs.py
    """

    root_directory = get_input_info()
    if root_directory == "None":
        print("Please provide a valid root directory path with --root_directory.")
    else:
        print(root_directory)

        for root, subdirectories, files in os.walk(root_directory):
            for file in files:
                if file == "config.yaml":
                    with open(os.path.join(root, file), 'r') as conf:
                        boolodejobs = bo.ConfigParser.parse(conf)
                    boolodejobs.execute_jobs()
    print('Done with running all config files.')

if __name__ == '__main__':
    """
     Calls the main function
    """
    main()