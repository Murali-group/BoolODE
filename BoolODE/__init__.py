import os
import sys
import yaml
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
from typing import Dict, List
# Local imports
from BoolODE import model_generator as mg
from BoolODE import run_experiment as runexp
from BoolODE import post_processing as po


class GlobalSettings(object):
    def __init__(self,
                 model_dir, output_dir,
                 do_simulations, do_post_processing,
                 modeltype) -> None:
        self.model_dir = model_dir
        self.output_dir = output_dir
        self.do_simulations = do_simulations
        self.do_post_processing = do_post_processing
        self.modeltype = modeltype

class JobSettings(object):
    '''
    TODO document
    '''
    def __init__(self, jobs) -> None:
        self.jobs = jobs
        
class PostProcSettings(object):
    def __init__(self,
                 dropout_jobs,
                 dimred_jobs,
                 slingshot_jobs) -> None:
        
        self.dropout_jobs = dropout_jobs
        self.dimred_jobs = dimred_jobs
        self.slingshot_jobs = slingshot_jobs


class BoolODE(object):
    '''
    The BoolODE object is created by parsing a user-provided configuration
    file. It then calls the correct functions for executing the user defined
    experiment.
    '''

    def __init__(self,
                 job_settings: JobSettings,
                 global_settings: GlobalSettings,
                 postproc_settings: PostProcSettings) -> None:
        self.job_settings = job_settings
        self.global_settings = global_settings
        self.post_settings = postproc_settings        
        self.jobs: Dict[int, Dict] = self.__process_jobs()

    def __process_jobs(self) -> Dict[int, Dict]:
        '''
        Creates a list of jobs, where each job is specified by 
        the values in a dictionary called data.
        Default parameter values are specified here.
        '''
        jobs = {}
        for jobid, job in enumerate(self.job_settings.jobs):
            data = {}
            # Create output folder if it doesnt exist
            data['name'] = job.get('name')
            data['outprefix'] = Path(self.global_settings.output_dir, job.get('name'))
            data['modelpath'] = Path(self.global_settings.model_dir, job.get('model_definition',''))
            data['simulation_time'] = job.get('simulation_time',20)
            data['icsPath'] = Path(self.global_settings.model_dir, job.get('model_initial_conditions',''))
            data['num_cells'] = job.get('num_cells',100)
            data['sample_cells'] = job.get('sample_cells',False)
            data['nClusters'] = job.get('nClusters',1)
            data['doParallel'] = job.get('do_parallel',False)            
            data['identical_pars'] = job.get('identical_pars',False)
            data['sample_pars'] = job.get('sample_pars',False)
            data['sample_std'] = job.get('sample_std',0.1)
            data['integration_step_size'] = job.get('integration_step_size',0.01)            
            # Optional Settings
            data['parameter_inputs_path'] = Path(self.global_settings.model_dir,\
                                                 job.get('parameter_inputs_path',''))
            data['parameter_set_path'] = Path(self.global_settings.model_dir,\
                                              job.get('parameter_set_path',''))
            data['interaction_strength_path'] = Path(self.global_settings.model_dir,\
                                                     job.get('interaction_strength_path',''))
            data['species_type_path'] = Path(self.global_settings.model_dir,\
                                             job.get('species_type_path',''))            
            # Simulator settings
            data['burnin'] = job.get('burn_in',False)
            data['writeProtein'] = job.get('write_protein',False)
            data['normalizeTrajectory'] = job.get('normalize_trajectory',False)
            data['add_dummy'] = job.get('add_dummy',False)
            data['max_parents'] = job.get('max_parents',1)
            data['modeltype'] = self.global_settings.modeltype

            jobs[jobid] = data
        return(jobs)

    def execute_jobs(self, parallel=False, num_threads=1):
        '''
        Run each user specified job. 
        BoolODE runs two types of functions
        1. If `do_simulation == TRUE`, perform SDE simulations of model specified as Boolean rules. 
        2. If `do_post_processing == TRUE` perform the list of post processing operations specified.

        Note:
        -----
        We recommend 
        
        Warning:
        -------
        This function automatically creates folders for each job name 
        as specified in the config file, if the folder doesn't already exist.
        Contents of existing folders will be rewritten!
        '''
        base_output_dir = self.global_settings.output_dir

        alljobs =  self.jobs.keys()
        print('Creating output folders')
        for jobid in alljobs:
            outdir = self.jobs[jobid]['outprefix']
            if not os.path.exists(outdir):
                print(outdir, "does not exist, creating it...")
                os.makedirs(outdir)
        if self.global_settings.do_simulations:
            print('Starting simulations')
            for jobid in alljobs:
                runexp.startRun(self.jobs[jobid])
        if self.global_settings.do_post_processing:
            print('Starting post processing')
            self.do_post_processing()

    def do_post_processing(self):
        """
        Call runSlingShot, and generateDropouts, if specified by the user.
        """
        alljobs =  self.jobs.keys()
        filetypedict = {'expr':'ExpressionData.csv',
                        'pseudo':'PseudoTime.csv',
                        'refNet':'refNetwork.csv'}        
        if self.post_settings.dropout_jobs is not None:
            print('Starting genDropouts...')
            for drop in self.post_settings.dropout_jobs:
                num_invalid = 0
                for jobid in alljobs:
                    settings = {}
                    invalid = False
                    settings['outPrefix'] = self.global_settings.output_dir +\
                        '/' + self.jobs[jobid]['name'] +\
                        '/' + self.jobs[jobid]['name']
                    settings['expr'] = Path(self.jobs[jobid]['outprefix']\
                                            ,'ExpressionData.csv')
                    settings['pseudo'] = Path(self.jobs[jobid]['outprefix']\
                                              ,'PseudoTime.csv')
                    settings['refNet'] = Path(self.jobs[jobid]['outprefix']\
                                              ,'refNetwork.csv')
                    settings['dropout'] = drop.get('dropout', True)                
                    settings['nCells'] = drop.get('nCells', 100)
                    settings['drop_cutoff'] = drop.get('drop_cutoff', 0.0)
                    settings['drop_prob'] = drop.get('drop_prob', 0.0)
                    
                    for filetype in ['expr', 'pseudo', 'refNet']:
                        if not settings[filetype].is_file():
                            print(self.jobs[jobid]['name'], ': ',filetypedict[filetype], "not found. Retry with `do_simulations: True` in global_settings.")                            
                            invalid = True
                            num_invalid += 1
                            break
                    if not invalid:
                        po.genDropouts(settings)
                if num_invalid == len(alljobs):
                    break
                    
        if self.post_settings.dimred_jobs is not None:
            print("Starting dimesionality reduction using tSNE")
            for dimred_jobs in self.post_settings.dimred_jobs:
                num_invalid = 0
                for jobid in alljobs:
                    settings = {}
                    invalid = False                    
                    settings['expr'] = Path(self.jobs[jobid]['outprefix'],\
                                            'ExpressionData.csv')
                    settings['pseudo'] = Path(self.jobs[jobid]['outprefix'],\
                                            'PseudoTime.csv')
                    settings['perplexity'] = dimred_jobs['perplexity']
                    for filetype in ['expr', 'pseudo']:
                        if not settings[filetype].is_file():
                            print(self.jobs[jobid]['name'], ':',filetypedict[filetype], "not found. Retry with `do_simulations: True` in global_settings.")                            
                            invalid = True
                            num_invalid += 1
                            break
                    if not invalid:
                        po.doDimRed(settings)
                        all_invalid = False
                if num_invalid == len(alljobs):
                    break
                        
        if self.post_settings.slingshot_jobs is not None:
            if self.post_settings.dimred_jobs is not None:
                print("Using default perplexity=300. (Specify `perplexity` under DimRed.)")
            print('Starting SlingShot...')
            for sshot in self.post_settings.slingshot_jobs:
                for jobid in alljobs:
                    settings = {}
                    invalid = False
                    settings['outPrefix'] = self.global_settings.output_dir +\
                        '/' + self.jobs[jobid]['name'] +\
                        '/' + self.jobs[jobid]['name'] + '-ss' 
                    settings['expr'] = Path(self.jobs[jobid]['outprefix']\
                                            ,'ExpressionData.csv')
                    settings['pseudo'] = Path(self.jobs[jobid]['outprefix']\
                                              ,'PseudoTime.csv')
                    settings['refNet'] = Path(self.jobs[jobid]['outprefix']\
                                              ,'refNetwork.csv')
                    if self.jobs[jobid]['nClusters'] == 1:
                        settings['nClusters'] = 1
                    else:
                        settings['nClusters'] = self.jobs[jobid]['nClusters'] + 1
                    settings['noEnd'] = sshot.get('noEnd', False)
                    settings['perplexity'] = sshot.get('perplexity', 300)

                    for filetype in ['expr', 'pseudo', 'refNet']:
                        if not settings[filetype].is_file():
                            print(self.jobs[jobid]['name'], ': ',filetypedict[filetype], "not found. Retry with `do_simulations: True` in global_settings.")
                            invalid = True
                            break
                    if not invalid:
                        po.computeSSPT(settings)                    
            
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
            ConfigParser.__parse_postproc_settings(
                config_map['post_processing'])            
        )

    @staticmethod
    def __parse_job_settings(joblist) -> JobSettings:
        jobs = joblist
        return JobSettings(jobs)
    
    @staticmethod
    def __parse_global_settings(input_settings_map) -> GlobalSettings:
        required_fields = ['model_dir',
                           'output_dir',
                           'do_simulations',
                           'do_post_processing',
                           'modeltype']
        print("Global settings:")
        print("----------------")        
        for rf in required_fields:
            try:
                print(f"\t{rf}: {input_settings_map[rf]}")
            except Exception:
                print("!!Missing inputs!!")
                print(f"BoolODE global settings is missing '{rf}'")
                print(f"Please specify {rf} in the input config file!\nExiting BoolODE...")
                sys.exit()
                
        model_dir = input_settings_map['model_dir']
        output_dir = input_settings_map['output_dir']
        do_simulations = input_settings_map['do_simulations']
        do_post_processing = input_settings_map['do_post_processing']
        modeltype = input_settings_map['modeltype']



        
        return GlobalSettings(model_dir,
                              output_dir,
                              do_simulations,
                              do_post_processing,
                              modeltype)
    @staticmethod
    def __parse_postproc_settings(input_settings_map) -> GlobalSettings:
        dropout_jobs = input_settings_map.get('Dropouts', None)
        slingshot_jobs = input_settings_map.get('Slingshot', None)
        dimred_jobs = input_settings_map.get('DimRed', None)
        return PostProcSettings(dropout_jobs, dimred_jobs, slingshot_jobs)        
