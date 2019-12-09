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
from BoolODE import simulatorUtils as  utils
from BoolODE import simulatorCore as simulator 
from importlib.machinery import SourceFileLoader
import multiprocessing as mp

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
                  interactionStrengthDF,
                  modeltype):
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
    heavisideSigma = 5

    # Chosen arbitrarily
    signalingtimescale = 5.0
    hillCoefficient = 10
    interactionStrength = 1.0

    parameterNamePrefixAndDefaultsAll = {
        # Hill coefficients
        'n_':hillCoefficient,
        # Thresholds
        'k_':hillThreshold,
        # Soft-Heaviside function steepness
        'sigmaH_':heavisideSigma
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
    ## TODO Add appropriate parameters if heaviside

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
                      interactionStrengthDF,
                      modeltype='hill'):
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
    :param modeltype: Specify type of modeling framework. ['hill', 'heaviside']
    :type modeltype: str
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

    par, x_max, parameterInputs = getParameters(DF, identicalPars,
                                                samplePars,
                                                sampleStd,
                                                genelist,
                                                proteinlist,
                                                withoutRules,
                                                parameterInputsDF,
                                                parameterSetDF,
                                                interactionStrengthDF,
                                                modeltype)

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
        if modeltype == 'hill':
            par['alpha_'+row['Gene']] = int(boolodespace['booleval'])
        elif modeltype == 'heaviside':
            par['omega_' + row['Gene']] = utils.heavisideThreshold(boolodespace['booleval'])

    for i,row in DF.iterrows():
        ## Parse Boolean rule to get list of regulators

        allreg, regSpecies, regInputs = utils.getRegulatorsInRule(row['Rule'],
                                                            species,
                                                            inputs)

        
        currSp = row['Gene']
        if modeltype == 'hill':
            num = '( alpha_' + currSp
            den = '( 1'
        elif modeltype == 'heaviside':
           exponent = '-sigmaH_' + currSp +'*( omega_' + currSp

        strengthSpecified = False
        regulatorsWithStrength = set()
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
            for combinationOfRegulators in combinations(allreg,i):
                regulatorExpression = utils.createRegulatoryTerms(combinationOfRegulators,
                                                                  modeltype,
                                                                  strengthSpecified,
                                                                  regulatorsWithStrength,
                                                                  regSpecies)

                if modeltype == 'hill':
                    # Create Numerator and Denominator
                    den += ' +' +  regulatorExpression
                    num += ' + a_' + currSp +'_'  + '_'.join(list(combinationOfRegulators)) + '*' + regulatorExpression
                elif modeltype == 'heaviside':
                    exponent += ' + w_' + currSp + '_' + '_'.join(list(combinationOfRegulators)) +'*' + regulatorExpression

                # evaluate rule to assign values to parameters
                ##################################################
                for i1, row1 in DF.iterrows():                 #
                    exec(row1['Gene'] + ' = 0', boolodespace)  #
                # Set each regulator to ON, evaluate rule. we are looping
                # over all such combinations of regulators
                for geneInList in combinationOfRegulators:     #
                    exec(geneInList + ' = 1', boolodespace)    #
                                                               #
                exec('boolval = ' + row['Rule'], boolodespace) #
                ##################################################

                if modeltype == 'hill':
                    par['a_' + currSp +'_'  + '_'.join(list(combinationOfRegulators))] = int(boolodespace['boolval'])
                elif modeltype == 'heaviside':
                    par['w_' + currSp +'_'  + '_'.join(list(combinationOfRegulators))] = 0.1*utils.heavisideThreshold(boolodespace['boolval'])                    

        if modeltype == 'hill':
            num += ' )'
            den += ' )'
            f = '(' + num + '/' + den + ')'
        elif modeltype == 'heaviside':
            exponent += ')'
            f = '(1./(1. + np.exp(np.sign('+exponent+')*min(100.,abs(' + exponent+ ')))))'
        
        if currSp in proteinlist:
            Production =  f
            Degradation = 'p_' + currSp
            varspecs['p_' + currSp] = 'signalingtimescale*(y_max*' + Production \
                                       + '-' + Degradation + ')'
        else:
            
            Production = 'm_'+ currSp + '*' + f
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
