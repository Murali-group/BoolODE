#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'
import os
import sys
import ast
import time
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from optparse import OptionParser
from itertools import combinations
from scipy.integrate import odeint
from sklearn.cluster import KMeans
from importlib.machinery import SourceFileLoader
import multiprocessing as mp
# local imports
from BoolODE import utils
from BoolODE import model_generator as mg
from BoolODE import simulator 

np.seterr(all='raise')

def Experiment(Model, ModelSpec,tspan,
               num_cells,
               sampleCells,
               varmapper, parmapper,
               genelist, proteinlist,
               outPrefix,icsDF,
               nClusters,
               x_max,
               doParallel,
               burnin=False,writeProtein=False,
               normalizeTrajectory=False):
    """
    Carry out an `in-silico` experiment. This function takes as input 
    an ODE model defined as a python function and carries out stochastic
    simulations. BoolODE defines a _cell_ as a single time point from 
    a simulated time course. Thus, in order to obtain 50 single cells,
    BoolODE carries out 50 simulations, which are stored in ./simulations/.
    Further, if it is known that the resulting single cell dataset will
    exhibit multiple trajectories, the user can specify  the number of clusters in
    `nClusters`; BoolODE will then cluster the entire simulation, such that each
    simulated trajectory possesess a cluster ID.

    :param Model: Function defining ODE model
    :type Model: function
    :param ModelSpec: Dictionary defining ODE model. See readBooleanRules()
    :type ModelSpec: dict
    :param tspan: Array of time points
    :type tspan: ndarray
    :param num_cells: Number of simulations to perform
    :type num_cells: int
    :param sampleCells: Bool that specifies if a random sample of size num_cells should be generated from the simulated output, where one cell is picked per simulation without replacement
    :type sampleCells: bool
    :param varmapper: 
    :type varmapper: dict
    :param parmapper: 
    :type parmapper: dict
    :param genelist: List of all gene names
    :type genelist:  list
    :param proteinlist: List of all protein names
    :type proteinlist: list
    :param outPrefix: Name of output folder. 
    :type outPrefix: str
    :param icsDF: Dataframe specifying initial condition for simulation
    :type icsDF: pandas DataFrame
    :param nClusters: Number of expected trajectories. Used to perform k-means clustering
    :type nClusters: int
    :param x_max: max value of gene. By default 2.0, will vary if parameters are sampled
    :type x_max: float
    :param doParallel: Bool specifying starting simulations in parallel
    :type doParallel: bool
    :param burnin: Bool specifying that initial fraction of simulation should be discarded. Obsolete
    :type burnin: bool
    :param writeProtein: Bool specifying if the protein values should be written to file. Default = False
    :type writeProtein: bool
    :param normalizeTrajectory: Bool specifying if the gene expression values should be scaled between 0 and 1.
    :type normalizeTrajectory: bool 
    """
    # if not sampleCells:
    #     print("Note: Simulated trajectories will be clustered. nClusters = %d" % nClusters)
    ####################    
    allParameters = dict(ModelSpec['pars'])
    parNames = sorted(list(allParameters.keys()))
    ## Use default parameters 
    pars = [ModelSpec['pars'][k] for k in parNames]
    ####################
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    revvarmapper = {v:k for k,v in varmapper.items()}
    proteinIndex = [i for i in range(len(varmapper.keys())) if 'p_' in varmapper[i]]

    y0 = [ModelSpec['ics'][varmapper[i]] for i in range(len(varmapper.keys()))]
    ss = np.zeros(len(varmapper.keys()))
    
    for i,k in varmapper.items():
        if 'x_' in k:
            ss[i] = 1.0
        elif 'p_' in k:
            if k.replace('p_','') in proteinlist:
                # Seting them to the threshold
                # causes them to drop to 0 rapidly
                # TODO: try setting to threshold < v < y_max
                ss[i] = 20.
            
    if icsDF is not None:
        icsspec = icsDF.loc[0]
        genes = ast.literal_eval(icsspec['Genes'])
        values = ast.literal_eval(icsspec['Values'])
        icsmap = {g:v for g,v in zip(genes,values)}
        for i,k in varmapper.items():
            for p in proteinlist:
                if p in icsmap.keys():
                    ss[revvarmapper['p_'+p]] = icsmap[p]
                else:
                    ss[revvarmapper['p_'+p]] = 0.01
            for g in genelist:
                if g in icsmap.keys():
                    ss[revvarmapper['x_'+g]] = icsmap[g]
                else:
                    ss[revvarmapper['x_'+g]] = 0.01
            
    if len(proteinlist) == 0:
        result = pd.DataFrame(index=pd.Index([varmapper[i] for i in rnaIndex]))
    else:
        speciesoi = [revvarmapper['p_' + p] for p in proteinlist]
        speciesoi.extend([revvarmapper['x_' + g] for g in genelist])
        result = pd.DataFrame(index=pd.Index([varmapper[i] for i in speciesoi]))
        
    # Index of every possible time point. Sample from this list
    startat = 0
    timeIndex = [i for i in range(startat, len(tspan))]        
    if sampleCells:
        # pre-define the time points from which a cell will be sampled
        # per simulation
        sampleAt = np.random.choice(timeIndex, size=num_cells)
        header = ['E' + str(cellid) + '_' + str(time) \
                  for cellid, time in\
                  zip(range(num_cells), sampleAt)]

    ## Construct dictionary of arguments to be passed
    ## to simulateAndSample(), done in parallel
    outPrefix = str(outPrefix)
    argdict = {}
    argdict['allParameters'] = allParameters
    argdict['parNames'] = parNames
    argdict['Model'] = Model
    argdict['tspan'] = tspan
    argdict['varmapper'] = varmapper
    argdict['timeIndex'] = timeIndex
    argdict['genelist'] = genelist
    argdict['proteinlist'] = proteinlist
    argdict['writeProtein'] = writeProtein
    argdict['outPrefix'] = outPrefix
    argdict['sampleCells'] = sampleCells
    argdict['pars'] = pars
    argdict['ss'] = ss
    argdict['ModelSpec'] = ModelSpec
    argdict['rnaIndex'] = rnaIndex
    argdict['proteinIndex'] = proteinIndex
    argdict['genelist'] = genelist
    argdict['proteinlist'] = proteinlist
    argdict['revvarmapper'] = revvarmapper
    argdict['x_max'] = x_max

    if sampleCells:
        argdict['header'] = header

    simfilepath = Path(outPrefix, './simulations/')
    if not os.path.exists(simfilepath):
        print(simfilepath, "does not exist, creating it...")
        os.makedirs(simfilepath)
    print('Starting simulations')
    start = time.time()
    ########################
    if doParallel:
        ### Contributed by matthieubulte
        with mp.Pool() as pool:
            jobs = []
            for cellid in range(num_cells):
                cell_args = dict(argdict, seed=cellid, cellid=cellid)
                job = pool.apply_async(simulateAndSample, args=(cell_args,))
                jobs.append(job)
                
            for job in jobs:
                job.wait()
        ###
    else:
        for cellid in tqdm(range(num_cells)):
            argdict['seed'] = cellid
            argdict['cellid'] = cellid
            simulateAndSample(argdict)
    ########################
    # ## Sleep for 1 s to allow for IO. Hack necessary for smalle number of simulations
    # ## where sim time < IO time.
    # time.sleep(1)
    print("Simulations took %0.3f s"%(time.time() - start))
    frames = []
    print('starting to concat files')
    start = time.time()
    if not sampleCells:
        # initialize dictionary to hold raveled values, used to cluster
        groupedDict = {} 
    for cellid in tqdm(range(num_cells)):
        if sampleCells:
            df = pd.read_csv(outPrefix + '/simulations/E'+str(cellid) + '-cell.csv',index_col=0)
            df = df.sort_index()                
        else:
            df = pd.read_csv(outPrefix + '/simulations/E'+str(cellid) + '.csv',index_col=0)
            df = df.sort_index()
            groupedDict['E' + str(cellid)] = df.values.ravel()
        frames.append(df.T)
    stop = time.time()
    print("Concating files took %.2f s" %(stop-start))
    result = pd.concat(frames,axis=0)
    result = result.T
    indices = result.index
    newindices = [i.replace('x_','') for i in indices]
    result.index = pd.Index(newindices)
    
    if nClusters > 1:
        ## Carry out k-means clustering to identify which
        ## trajectory a simulation belongs to
        print('Starting k-means clustering')
        groupedDF = pd.DataFrame.from_dict(groupedDict)
        print('Clustering simulations...')
        start = time.time()            
        # Find clusters in the experiments
        clusterLabels= KMeans(n_clusters=nClusters,
                              n_jobs=8).fit(groupedDF.T.values).labels_
        print('Clustering took %0.3fs' % (time.time() - start))
        clusterDF = pd.DataFrame(data=clusterLabels, index =\
                                 groupedDF.columns, columns=['cl'])
        clusterDF.to_csv(outPrefix + '/ClusterIds.csv')
    else:
        print('Requested nClusters=1, not performing k-means clustering')
    ##################################################
    
    return result
    
def startRun(settings):
    validInput = utils.checkValidModelDefinitionPath(settings['path'], settings['name'])
    if not validInput:
        return
    
    startfull = time.time()

    outdir = settings['outprefix']
    if not os.path.exists(outdir):
        print(outdir, "does not exist, creating it...")
        os.makedirs(outdir)

    parameterInputsDF = utils.checkValidInputPath(settings['parameter_inputs_path'])
    parameterSetDF = utils.checkValidInputPath(settings['parameter_set_path'])
    icsDF = utils.checkValidInputPath(settings['icsPath'])
    interactionStrengthDF = utils.checkValidInputPath(settings['interaction_strength_path']) 
    speciesTypeDF = utils.checkValidInputPath(settings['species_type_path'])
    
    tmax = settings['simulation_time']    
    integration_step_size = settings['integration_step_size']
    tspan = np.linspace(0,tmax,int(tmax/integration_step_size))

    modelgenerator = ModelGenerator(settings)
    # rulesdf, withoutRules = mg.readBooleanRules(settings['path'],
    #                                          settings['parameter_inputs_path'],
    #                                          settings['outprefix'],
    #                                          settings['add_dummy'],
    #                                          settings['max_parents'])

    it = 0
    someexception = True
    while someexception:
        try:
            genesDict = {}
            
            ModelSpec,\
                parameterInputs,\
                genelist,\
                proteinlist,\
                x_max = mg.generateModelDict(rulesdf,
                                          settings['identical_pars'],
                                          settings['sample_pars'],
                                          settings['sample_std'],
                                          withoutRules,
                                          speciesTypeDF,
                                          parameterInputsDF,
                                          parameterSetDF,
                                          interactionStrengthDF,
                                          settings['modeltype'])
            # FIXME : ask user to pass row if parameter input file
            # Hardcoded. We only care about one input vector, say the first one
            # TODO check this
            # this seems unnecessary
            if parameterInputsDF is not None:
                ModelSpec['pars'].update(parameterInputs[0])
                
            varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
            parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
            path_to_model  = utils.writeModelToFile(ModelSpec, outPrefix=settings['outprefix'])
            utils.writeParametersToFile(ModelSpec, settings['outprefix'])            
            ## Load file from path
            print(path_to_model.as_posix())
            model = SourceFileLoader("model", path_to_model.as_posix()).load_module()
            resultDF = Experiment(model.Model,
                                  ModelSpec,
                                  tspan,
                                  settings['num_cells'],
                                  settings['sample_cells'],
                                  varmapper, parmapper,
                                  genelist, proteinlist,
                                  settings['outprefix'],
                                  icsDF,
                                  settings['nClusters'],
                                  x_max,
                                  settings['doParallel'],
                                  burnin=settings['burnin'],
                                  writeProtein=settings['writeProtein'],
                                  normalizeTrajectory=settings['normalizeTrajectory'])
            print('Generating input files for pipline...')
            start = time.time()
            utils.generateInputFiles(resultDF, rulesdf,
                                     withoutRules,
                                     parameterInputsDF,tmax,settings['num_cells'],
                                     outPrefix=settings['outprefix'])
            print('Input file generation took %0.2f s' % (time.time() - start))
            someexception= False
            
        except FloatingPointError as e:
            it +=1 
            print(e,"\nattempt %d" %it)
        

    print("BoolODE.py took %0.2fs"% (time.time() - startfull))
    print('all done.')    

def simulateAndSample(argdict):
    allParameters = argdict['allParameters']
    parNames = argdict['parNames']
    Model = argdict['Model']
    tspan = argdict['tspan']
    varmapper = argdict['varmapper']
    timeIndex = argdict['timeIndex']
    genelist = argdict['genelist']
    proteinlist = argdict['proteinlist']
    writeProtein=argdict['writeProtein']
    cellid = argdict['cellid']
    outPrefix = argdict['outPrefix']
    sampleCells = argdict['sampleCells']
    ss = argdict['ss']
    ModelSpec = argdict['ModelSpec']
    rnaIndex = argdict['rnaIndex']
    proteinIndex = argdict['proteinIndex']
    genelist = argdict['genelist']
    proteinlist = argdict['proteinlist']
    revvarmapper = argdict['revvarmapper']
    seed = argdict['seed']
    pars = argdict['pars']
    x_max = argdict['x_max']
    
    # Retained for debugging
    isStochastic = True
    
    if sampleCells:
        header = argdict['header']
        
    pars = {}
    for k, v in allParameters.items():
        pars[k] = v
    pars = [pars[k] for k in parNames]
    
    ## Boolean to check if a simulation is going to a
    ## 0 steady state, with all genes/proteins dying out
    retry = True
    trys = 0
    ## timepoints
    tps = [i for i in range(1,len(tspan))]
    ## gene ids
    gid = [i for i,n in varmapper.items() if 'x_' in n]
    outPrefix = outPrefix + '/simulations/'
    while retry:
        seed += 1000
        y0_exp = mg.getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex,
                                     genelist, proteinlist,
                                     varmapper,revvarmapper)
        
        P = simulator.simulateModel(Model, y0_exp, pars, isStochastic, tspan, seed)
        P = P.T
        retry = False
        ## Extract Time points
        subset = P[gid,:][:,tps]
        df = pd.DataFrame(subset,
                          index=pd.Index(genelist),
                          columns = ['E' + str(cellid) +'_' +str(i)\
                                     for i in tps])
        ## Heuristic:
        ## If the largest value of a protein achieved in a simulation is
        ## less than 10% of the y_max, drop the simulation.
        ## This check stems from the observation that in some simulations,
        ## all genes go to the 0 steady state in some rare simulations.
        dfmax = df.max()
        for col in df.columns:
            colmax = df[col].max()
            if colmax < 0.1*x_max:
                retry= True
                break
        
        if sampleCells:
            ## Write a single cell to file
            ## These samples allow for quickly and
            ## reproducibly testing the output.
            sampledf = utils.sampleCellFromTraj(cellid,
                                          tspan, 
                                          P,
                                          varmapper, timeIndex,
                                          genelist, proteinlist,
                                          header,
                                          writeProtein=writeProtein)
            sampledf = sampledf.T
            sampledf.to_csv(outPrefix + 'E' + str(cellid) + '-cell.csv')            
            
        trys += 1
        # write to file
        df.to_csv(outPrefix + 'E' + str(cellid) + '.csv')
        
        if trys > 1:
            print('try', trys)
            

    
