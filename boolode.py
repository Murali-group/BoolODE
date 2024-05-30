import argparse
import BoolODE as bo

def get_parser():
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Run BoolODE to generate synthetic scRNAseq data.\n Please specify config file.')

    parser.add_argument('--config', default='config.yaml',
        help='Path to config file')

    return parser

def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts

def main():
    opts = parse_arguments()
    config_file = opts.config
    with open(config_file, 'r') as conf:
        boolodejobs = bo.ConfigParser.parse(conf)
    
    boolodejobs.execute_jobs()
    print('Jobs finished')


if __name__ == '__main__':
    main()

