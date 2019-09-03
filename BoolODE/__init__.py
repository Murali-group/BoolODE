import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
from typing import Dict, List
from BoolODE import model_generator as mg
import os

class GlobalSettings(object):
    def __init__(self,
                 model_dir, output_dir) -> None:
        
        self.model_dir = model_dir
        self.output_dir = output_dir

class JobSettings(object):
    '''
    TODO document
    '''

    def __init__(self, jobs) -> None:
        self.jobs = jobs

class BoolODE(object):
    '''
    The BoolODE object is created by parsing a user-provided configuration
    file. It then calls the correct functions for executing the user defined
    experiment.
    '''

    def __init__(self,
            job_settings: JobSettings,
            global_settings: GlobalSettings) -> None:

        self.job_settings = job_settings
        self.global_settings = global_settings
        self.jobs: Dict[int, Dict] = self.__process_jobs()

    def __process_jobs(self) -> Dict[int, Dict]:
        '''
        TODO add docs
        '''
        jobs = {}#: Dict[int, Dict] = defaultdict(list)
        for jobid, job in enumerate(self.job_settings.jobs):
            data = {}
            # Create output folder if it doesnt exust
            data['name'] = job.get('name')
            data['outprefix'] = Path(self.global_settings.output_dir, job.get('name'))
            data['path'] = Path(self.global_settings.model_dir, job.get('model_definition',''))
            data['simulation_time'] = job.get('simulation_time',20)
            data['icsPath'] = Path(self.global_settings.model_dir, job.get('model_initial_conditions',''))
            data['num_cells'] = job.get('num_cells',100)
            data['sample_cells'] = job.get('sample_cells',False)
            data['nClusters'] = job.get('nClusters',1)
            data['doParallel'] = job.get('do_parallel',False)            
            data['identical_pars'] = job.get('identical_pars',False)
            data['sample_pars'] = job.get('sample_pars',False)
            data['sample_std'] = job.get('std',0.1)
            data['integration_step_size'] = job.get('integration_step_size',0.01)            
            # Optional Settings
            data['parameter_inputs_path'] = Path(self.global_settings.model_dir, job.get('parameter_inputs_path',''))
            data['parameter_set_path'] = Path(self.global_settings.model_dir, job.get('parameter_set_path',''))
            data['interaction_strength_path'] = Path(self.global_settings.model_dir, job.get('interaction_strength_path',''))
            data['species_type_path'] = Path(self.global_settings.model_dir, job.get('species_type_path',''))            
            # Simulator settings
            data['burnin'] = job.get('burn_in',False)
            data['writeProtein'] = job.get('write_protein',False)
            data['normalizeTrajectory'] = job.get('normalize_trajectory',False)

            data['add_dummy'] = job.get('add_dummy',False)
            data['max_parents'] = job.get('max_parents',1)

            jobs[jobid] = data
        return(jobs)

    def execute_jobs(self, parallel=False, num_threads=1):
        '''
        Run each of the algorithms
        '''
        base_output_dir = self.global_settings.output_dir

        alljobs =  self.jobs.keys()
        for jobid in alljobs:
            mg.startRun(self.jobs[jobid])
        
class ConfigParser(object):
    '''
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline
    '''
    @staticmethod
    def parse(config_file_handle) -> BoolODE:
        config_map = yaml.safe_load(config_file_handle)
        return BoolODE(
            ConfigParser.__parse_job_settings(
                config_map['jobs']),
            ConfigParser.__parse_global_settings(
                config_map['global_settings']),            
        )

    @staticmethod
    def __parse_job_settings(joblist) -> JobSettings:
        jobs = joblist
        
        return JobSettings(jobs)
    
    @staticmethod
    def __parse_global_settings(input_settings_map) -> GlobalSettings:
        model_dir = input_settings_map['model_dir']
        output_dir = input_settings_map['output_dir']

        return GlobalSettings(model_dir,
                              output_dir)    
    
    # @staticmethod
    # def __parse_input_settings(input_settings_map) -> InputSettings:
    #     model_dir = input_settings_map['model_dir']
    #     dataset_dir = input_settings_map['dataset_dir']
    #     datasets = input_settings_map['datasets']

    #     return InputSettings(
    #             Path(input_dir, dataset_dir),
    #             datasets,
    #             ConfigParser.__parse_algorithms(
    #             input_settings_map['algorithms']))


    # @staticmethod
    # def __parse_algorithms(algorithms_list):
    #     algorithms = []
    #     for algorithm in algorithms_list:
    #             combos = [dict(zip(algorithm['params'], val))
    #                 for val in itertools.product(
    #                     *(algorithm['params'][param]
    #                         for param in algorithm['params']))]
    #             for combo in combos:
    #                 algorithms.append([algorithm['name'],combo])
            

    #     return algorithms

    # @staticmethod
    # def __parse_output_settings(output_settings_map) -> OutputSettings:
    #     output_dir = Path(output_settings_map['output_dir'])
    #     output_prefix = Path(output_settings_map['output_prefix'])

    #     return OutputSettings(output_dir,
    #                          output_prefix)
