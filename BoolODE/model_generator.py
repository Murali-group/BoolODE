#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'
import os
import sys
import ast
import time
import importlib
import warnings
import numpy as np
import scipy as sc
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from optparse import OptionParser
from itertools import combinations
from scipy.integrate import odeint
from sklearn.cluster import KMeans
from BoolODE import simulatorUtils as  utils
from BoolODE import simulatorCore as simulator 
from importlib.machinery import SourceFileLoader
from multiprocessing import Process, Lock, Manager

np.seterr(all='raise')

def readBooleanRules(path, parameterInputsPath, outPrefix='',
                     add_dummy=False, max_parents=1):
    """
    Reads a rule file from path and stores the rules in a dictionary
    
    Parameters
    ----------
    :param path: Path to Boolean Rule file
    :type path: str
    :param parameterInputsPath: Path to file containing input parameters to the model
    :type parameterInputsPath: str
    :param add_dummy: Experimental feature, adds dummy targets to each gene in the model. Useful to change the density of the network
    :type add_dummy: bool
    :param max_parents: Experimental feature, specifies number of parent genes per dummy gene added
    :type max_parents: int

    :returns:
        - df: Dataframe containing rules
        - withoutrules: list of nodes in input file without rules
    """
    df = pd.read_csv(path,sep='\t',engine='python')
    withRules = list(df['Gene'].values)
    allnodes = set()
    for ind,row in df.iterrows():
        rhs = row['Rule']
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')

        reg = [t for t in tokens if t not in ['not','and','or','']]
        allnodes.update(set(reg))

    withoutRules = list(allnodes.difference(set(withRules)))
    for n in withoutRules:
        if len(parameterInputsPath) == 0:
            print(n, "has no rule, adding self-activation.")
            df = df.append({'Gene':n,'Rule':n},ignore_index=True)
        else:
            print("Treating %s as parameter" % {n})
    ## Add dummy genes
    ## determinisitcally add parents
    if add_dummy:
        print('using max_parents=' + str(max_parents))
        num_dummy = 2*len(allnodes)
        for dg in range(num_dummy):
            ## Variable number of parents
            # DF = DF.append({'Gene':'dummy' + str(dg),
            #                 'Rule':' or '.join([s for s in \
            #                                     np.random.choice(list(allnodes),
            #                                                      size=np.random.choice([i+1 for i in range(len(allnodes)-1)]))])},
            #                ignore_index = True)
            df = df.append({'Gene':'dummy' + str(dg),
                            'Rule':' or '.join([s for s in \
                                                np.random.choice(list(allnodes),
                                                                 size=max_parents)])},
                           ignore_index = True)
            df.to_csv(outPrefix + 'rules-with-added-genes.csv')
            
    print(df)
    return df, withoutRules

def getParameters(DF,identicalPars,
                  samplePars,
                  sampleStd,
                  genelist,
                  proteinlist,
                  withoutRules,
                  parameterInputsDF,
                  parameterSetDF,
                  interactionStrengthDF):
    """
    Create dictionary of parameters and values. Assigns
    parameter values by evaluating the Boolean expression
    for each variable.

    :param DF: Table of values with two columns, 'Gene' specifies target, 'Rule' specifies Boolean function
    :type DF: pandas DataFrame
    :param identicalPars: Passed to utils.getSaneNval to set identical parameters
    :type identicalPars: bool
    :param samplePars: Sample kinetic parameters using a Gaussian distribution  centered around the default parameters
    :type samplePars: bool
    :param parameterInputsDF: Optional table that specifies input parameter values. Useful to specify experimental conditions.
    :type parameterInputsDF: pandas DataFrame
    :param parameterSetDF: Optional table that specifies predefined parameter set.
    :type parameterSetDF: pandas DataFrame
    :param interactionStrengthDF: Optional table that specifies interaction strengths. When not specified, default strength is set to 1.
    :type interactionStrengthDF: pandas DataFrame
    :returns:
        pars : dict
            Dictionary of parameters 
    """
    ## Set default parameters
    mRNATranscription = 20.
    mRNADegradation = mRNATranscription/2.
    proteinTranslation = 10
    proteinDegradation = 1
    ## The threshold is calculated as the max value of the species
    ## divided by 2.
    y_max = (proteinTranslation/proteinDegradation)*\
            (mRNATranscription/mRNADegradation)
    x_max = (mRNATranscription/mRNADegradation)
    hillThreshold = y_max/2

    # Chosen arbitrarily
    signalingtimescale = 5.0
    hillCoefficient = 10
    interactionStrength = 1.0

    parameterNamePrefixAndDefaultsAll = {
        # Hill coefficients
        'n_':hillCoefficient,
        # Thresholds
        'k_':hillThreshold,
    }     
    parameterNamePrefixAndDefaultsGenes = {
        # mRNA transcription rates
        'm_':mRNATranscription,
        # mRNA degradation rates
        'l_x_':mRNADegradation,
        # Protein translation rate
        'r_':proteinTranslation,
        # protein degradation rates
        'l_p_':proteinDegradation
    }     
    
    par = dict()
    ## If there is even one signaling protein,
    ## create the y_max parameter
    if len(proteinlist) > 0:
        par['y_max'] = y_max
        par['signalingtimescale'] = signalingtimescale
        
    parameterInputsDict  = {}
    interactionStrengths = {}

    ## Get set of species which have rules specified for them
    species = set(DF['Gene'].values)
    
    ## Check 1:
    ## If a parameter input file is specified,
    ## Every row in parameterInputsDF contains two
    ## columns:
    ##     'Inputs' is a python list of strings
    ## specifying the parameters.
    ##     'Values' is a python list of floats
    ## each corresponding to the value of a parameter.
    ## Now, we create a fake hill term for this input, even
    ## though there is no equation corresponding to it. Thus,
    ## we create a hillThreshold and hillCoefficient 
    ## term for it. This is useful so we can treat this input
    ## parameter as a "regulator" and leave the structure of
    ## the logic term unchanged.
    if parameterInputsDF is not None:
        for i, row in parameterInputsDF.iterrows():
            # initialize
            parameterInputsDict[i] = {p:0 for p in withoutRules}
            inputParams = ast.literal_eval(row['Inputs'])
            inputValues = ast.literal_eval(row['Values'])
            # Set to a max value
            parameterInputsDict[i].update({p:hillThreshold*2 for p,v in\
                                       zip(inputParams,inputValues) if v > 0})
            
        # Add input parameters, set to 0 by default
        par.update({p:0 for p in parameterInputsDict[0].keys()})
        par.update({'k_'+p:hillThreshold for p in parameterInputsDict[0].keys()})
        par.update({'n_'+p:hillCoefficient for p in parameterInputsDict[0].keys()})


    ## Check 2:
    ## If interaction strengths are specified.
    ## Create a new parameter specifying the strength
    ## of a given interaction. For all other equations in
    ## which "regulator" appears in, its hillThreshold value
    ## will remain the default value.
    if interactionStrengthDF is not None:
        for i,row in interactionStrengthDF.iterrows():
            regulator = row['Gene2']
            target = row['Gene1']
            interactionStrength = row['Strength']
            par.update({'k_' + regulator + '_' + target: hillThreshold/interactionStrength})

    ## Check 3:
    ## Whether to sample parameters or stick to defaults
    if samplePars:
        print("Sampling parameters")
        sigmamult = sampleStd
        print("Using std=" + str(sampleStd))
        lomult = 0.9
        himult = 1.1
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            sampledParameterValues = utils.getSaneNval(len(species),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=sigmamult*parDefault,\
                                                 identicalPars=identicalPars)
            for sp, sparval in zip(species, sampledParameterValues):
                if sp in genelist:
                    par[parPrefix + sp] = sparval
                    
        transcriptionRate = 0.0
        mRNADegradationRate = 0.0
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            sampledParameterValues = utils.getSaneNval(len(species),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=sigmamult*parDefault,\
                                                 identicalPars=identicalPars)
            if identicalPars:
                if parPrefix == 'm_':
                    transcriptionRate = sampledParameterValues[0]
                if parPrefix == 'l_x_':
                    mRNADegradationRate = sampledParameterValues[0]
            else:
                if parPrefix == 'm_':
                    transcriptionRate = parDefault
                if parPrefix == 'l_x_':
                    mRNADegradationRate = parDefault
            
            for sp, sparval in zip(species, sampledParameterValues):
                if sp in genelist:
                    par[parPrefix + sp] = sparval
        ## Based on these samples, estimate max value a
        ## gene can take at _steady-state_
        x_max = (transcriptionRate/mRNADegradationRate)
    else:
        print("Fixing rate parameters to defaults")
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            for sp in species:
                par[parPrefix + sp] = parDefault

        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            for sp in species:
                if sp in genelist:
                    par[parPrefix + sp] = parDefault
    ## Final Check.
    ## If parameterSetDF is specified, reassign
    ## parameter values to those from table
    ## This guarantees that all the parameters are defined without
    ## specific checks in parameterSetDF
    if parameterSetDF is not None:
        for pname, pvalue in parameterSetDF.iterrows():
            par[pname] = pvalue[1]
    return par, x_max, parameterInputsDict

def generateModelDict(DF,identicalPars,
                      samplePars,
                      sampleStd,
                      withoutRules,
                      speciesTypeDF,
                      parameterInputsDF,
                      parameterSetDF,
                      interactionStrengthDF):
    """
    Take a DataFrame object with Boolean rules,
    construct ODE equations for each variable.
    Takes optional parameter inputs and interaction
    strengths

    :param DF: Table of values with two columns, 'Gene' specifies target, 'Rule' specifies Boolean function
    :type DF: pandas DataFrame
    :param identicalPars: Passed to utils.getSaneNval to set identical parameters
    :type identicalPars: bool
    :param samplePars: Sample kinetic parameters using a Gaussian distribution  centered around the default parameters
    :type samplePars: bool
    :param sampleStd: Sample from a distribution of mean=par, sigma=sampleStd*par
    :type sampleStd: float
    :param speciesTypeDF: Table defining type of each species. Takes values  'gene' or 'protein'
    :type speciesTypeDF: pandas DataFrame
    :param parameterInputsDF: Optional table that specifies parameter input values. Useful to specify experimental conditions.
    :type parameterInputsDF: pandas DataFrame
    :param parameterSetDF: Optional table that specifies a predefined parameter set.
    :type parameterSetDF: pandas DataFrame
    :param interactionStrengthDF: Optional table that specifies interaction strengths. When not specified, default strength is set to 1.
    :type interactionStrengthDF: pandas DataFrame
    :returns:
        ModelSpec : dict
            Dictionary of dictionaries. 
            'varspecs' {variable:ODE equation}
            'ics' {variable:initial condition}
            'pars' {parameter:parameter value}
    """
    species = set(DF['Gene'].values)

    # Variables:
    ## Check:
    ## If speciesType is provided, initialize the
    ## correct species names
    ## Rule:
    ## Every gene (x_species) gets a corresponding
    ## protein (p_species). The gene equation contains
    ## regulatory logic, protein equation contains production
    ## degradation terms.
    ## Every 'protein' (p_species) only gets one variable,
    ## and the corresponding equation contains the regulatory
    ## logic term. There is no production/degradation associated
    ## with this species.
    varspecs = {}
    genelist = []
    proteinlist = []
    if speciesTypeDF is None:
        # Assume everything is a gene.
        for sp in species:
            varspecs['x_' + sp] = ''
            varspecs['p_' + sp] = ''
            genelist.append(sp)
    else:
        for i, row in speciesTypeDF.iterrows():
            sp = row['Species']
            if row['Type'] == 'protein':
                varspecs['p_' + sp] = ''
                proteinlist.append(sp)
            elif row['Type'] == 'gene':
                varspecs['x_' + sp] = ''
                varspecs['p_' + sp] = ''
                genelist.append(sp)
        specified = set(speciesTypeDF['Species'].values)
        for sp in species:
            if sp not in specified:
                genelist.append(sp)
                
    ## Useful for debugging. 
    # print(genelist)
    # print(proteinlist)

    par, x_max, parameterInputs = getParameters(DF,identicalPars,
                                                samplePars,
                                                sampleStd,
                                                genelist,
                                                proteinlist,
                                                withoutRules,
                                                parameterInputsDF,
                                                parameterSetDF,
                                                interactionStrengthDF)

    ##########################################################
    ## Assign values to alpha parameters representing the logic
    ## relationships between variables
    # Initialize new namespace
    boolodespace = {}
    for i,row in DF.iterrows():
        # Initialize species to 0
        tempStr = row['Gene'] + " = 0"  
        exec(tempStr, boolodespace)

    if parameterInputsDF is None:
        inputs = set()
    else:
        inputs = set(withoutRules)
        for k in parameterInputs[0].keys():
            # Initialize variables to 0
            tempStr = k + " = 0"  
            exec(tempStr, boolodespace)
            
    for i,row in DF.iterrows():
        # Basal alpha:
        # Execute the rule to figure out
        # the value of alpha_0
        exec('booleval = ' + row['Rule'], boolodespace) 
        par['alpha_'+row['Gene']] = int(boolodespace['booleval'])

    for i,row in DF.iterrows():
        rhs = row['Rule']
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')

        allreg = set([t for t in tokens if (t in species or t in inputs)])
        regulatorySpecies = set([t for t in tokens if t in species])
        inputreg = set([t for t in tokens if t in inputs])
        
        currSp = row['Gene']
        num = '( alpha_' + currSp
        den = '( 1'

        strengthSpecified = False
        
        if interactionStrengthDF is not None:
            if currSp in interactionStrengthDF['Gene1'].values:
                strengthSpecified = True
                # Get the list of regulators of currSp
                # whose interaction strengths have been specified
                regulatorsWithStrength = set(interactionStrengthDF[\
                                                                   interactionStrengthDF['Gene1']\
                                                                   == currSp]['Gene2'].values)
                print(regulatorsWithStrength)
                
        for i in range(1,len(allreg) + 1):
            for c in combinations(allreg,i):
                # Create the hill function terms for each regulator
                hills = []
                for ci in c:
                    if strengthSpecified and ci in regulatorsWithStrength:
                            hillThresholdName = 'k_' + ci + '_' + currSp
                    else:
                        hillThresholdName = 'k_' + ci
                        
                    if ci in regulatorySpecies:
                        # Note: Only proteins can be regulators
                        hills.append('(p_'+ci+'/'+hillThresholdName+')^n_'+ci)
                    elif ci in inputreg:
                        hills.append('('+ci+'/'+hillThresholdName+')^n_'+ci)
                mult = '*'.join(hills)
                # Create Numerator and Denominator
                den += ' +' +  mult                
                num += ' + a_' + currSp +'_'  + '_'.join(list(c)) + '*' + mult
                
                for i1, row1 in DF.iterrows():
                    exec(row1['Gene'] + ' = 0', boolodespace)

                for geneInList in c:
                    exec(geneInList + ' = 1', boolodespace)
                exec('boolval = ' + row['Rule'], boolodespace)

                par['a_' + currSp +'_'  + '_'.join(list(c))] = int(boolodespace['boolval']) 

        num += ' )'
        den += ' )'
        
        if currSp in proteinlist:
            Production =  '(' +num + '/' + den + ')'
            Degradation = 'p_' + currSp
            varspecs['p_' + currSp] = 'signalingtimescale*(y_max*' + Production \
                                       + '-' + Degradation + ')'
        else:
            Production = 'm_'+ currSp + '*(' +num + '/' + den + ')'
            Degradation = 'l_x_'  + currSp + '*x_' + currSp
            varspecs['x_' + currSp] =  Production \
                                       + '-' + Degradation
            # Create the corresponding translated protein equation
            varspecs['p_' + currSp] = 'r_'+currSp+'*'+'x_' +currSp + '- l_p_'+currSp+'*'+'p_' + currSp
            
    ##########################################################                                       
        
    # Initialize variables between 0 and 1, Doesn't matter.
    xvals = [1. for _ in range(len(genelist))]
    pvals = [20. for _ in range(len(proteinlist))]    
    ics = {}

    for sp, xv in zip(genelist, xvals):
        ics['x_' + sp] = xv
        ics['p_' + sp] = 0
    for sp, pv in zip(proteinlist, pvals):
        ics['p_' + sp] = pv

    ModelSpec = {}
    ModelSpec['varspecs'] = varspecs
    ModelSpec['pars'] = par
    ModelSpec['ics'] = ics
    return ModelSpec, parameterInputs, genelist, proteinlist, x_max

def getInitialCondition(ss, ModelSpec, rnaIndex,
                        proteinIndex,
                        genelist, proteinlist,
                        varmapper,revvarmapper):
    """
    Calculate the initial values of all state variables. 
    Takes into consideration user defined initial conditions, and computes the steady 
    states of the protein variables based on the estimated values of their corresponding genes.

    :param ss: Steady state array
    :type ss: ndarray
    :param ModelSpec: Dictionary of dictionary specifying the ODE model, containing parameters, initial conditions and equations.
    :type ModelSpec: dict
    :param rnaIndex: list of indices of genes
    :type rnaIndex: list
    :param proteinIndex: List of indices of proteins
    :type proteinIndex: list
    :param genelist: List of names of all genes in the model
    :type genelist: list
    :param proteinlist: List of names of all proteins in the model
    :type proteinlist: list
    :param varmapper: Mapper: {variable name : index}
    :type varmapper: dict
    :param revvarmapper: Mapper: {index : variable name}
    :type revvarmapper: dict
    :returns:
        - newics: List containing new initial conditions
    """
    
    # Initialize
    new_ics = [0 for _ in range(len(varmapper.keys()))]
    # Set the mRNA ics
    for ind in rnaIndex:
        if ss[ind] < 0:
            ss[ind] = 0.0
        new_ics[ind] =  ss[ind]
        if new_ics[ind] < 0:
            new_ics[ind] = 0
    for p in proteinlist:
        ind = revvarmapper['p_'+p]
        if ss[ind] < 0:
            ss[ind] = 0.0
        new_ics[ind] =  ss[ind]
        if new_ics[ind] < 0:
            new_ics[ind] = 0            
    # Calculate the Protein ics based on mRNA levels
    for genename in genelist:
        pss = ((ModelSpec['pars']['r_' + genename])/\
                                      (ModelSpec['pars']['l_p_' + genename]))\
                                      *new_ics[revvarmapper['x_' + genename]]
        new_ics[revvarmapper['p_' + genename.replace('_','')]] = pss
    return(new_ics)

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
        y0_exp = getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex,
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
        if trys > 1:
            print('try', trys)
            

    # write to file
    df.to_csv(outPrefix + 'E'+ str(cellid) + '.csv')
    
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
    if not sampleCells:
        print("Note: Simulated trajectories will be clustered. nClusters = %d" % nClusters)
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
        ## Carry out simulations in parrallel
        lock = Lock()
        jobs = []
        for cellid in range(num_cells):
            argdict['seed'] = cellid
            argdict['cellid'] = cellid
            p = Process(target=simulateAndSample,args=([argdict]))
            p.start()
            #p.join()
    else:
        for cellid in tqdm(range(num_cells)):
            argdict['seed'] = cellid
            argdict['cellid'] = cellid
            simulateAndSample(argdict)
    ########################
    ## Sleep for 1 s to allow for IO. Hack necessary for smalle number of simulations
    ## where sim time < IO time.
    time.sleep(1)
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
    
    if not sampleCells:
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
    
    rulesdf, withoutRules = readBooleanRules(settings['path'],
                                             settings['parameter_inputs_path'],
                                             settings['outprefix'],
                                             settings['add_dummy'],
                                             settings['max_parents'])

    # TODO cleanup
    # if len(withoutRules) == 0:
    #     withoutRules = []
        
    it = 0
    someexception = True
    while someexception:
        try:
            genesDict = {}
            
            ModelSpec,\
                parameterInputs,\
                genelist,\
                proteinlist,\
                x_max = generateModelDict(rulesdf,
                                          settings['identical_pars'],
                                          settings['sample_pars'],
                                          settings['sample_std'],
                                          withoutRules,
                                          speciesTypeDF,
                                          parameterInputsDF,
                                          parameterSetDF,
                                          interactionStrengthDF)
            # FIXME : ask user to pass row if parameter input file
            # Hardcoded. We only care about one input vector, say the first one
            # TODO check this
            # this seems unnecessary
            if parameterInputsDF is not None:
                ModelSpec['pars'].update(parameterInputs[0])
                
            varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
            parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
            dir_path  = utils.writeModelToFile(ModelSpec)
            ## Load file from path
            model = SourceFileLoader("model", dir_path + "/model.py").load_module()
            #import model
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
                                     parameterInputsDF,
                                     outPrefix=settings['outprefix'])
            print('Input file generation took %0.2f s' % (time.time() - start))
            
            someexception= False
            
        except FloatingPointError as e:
            it +=1 
            print(e,"\nattempt %d" %it)
        
    utils.writeParametersToFile(ModelSpec, settings['outprefix'])
    print("BoolODE.py took %0.2fs"% (time.time() - startfull))
    print('all done.')    

