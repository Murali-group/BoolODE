import os
import sys
import yaml
import numpy as np
import pandas as pd
from pathlib import Path


def heavisideThreshold(value):
    """
    Decides whether the basal parameter omega_0 should 
    be 1 (basal expression) or -1 (basal repression)
    """
    if int(value) == 1:
        return 1
    elif int(value) == 0:
        return -1

def loadParameterValues():
    """
    Checks for valid parameters in parameters.yaml file
    """
    with open(str(Path(__file__).parent.absolute()) + '/parameters.yaml','r') as parameterfile:
        parameterDefaults = yaml.safe_load(parameterfile)
        
    requiredPars = ['mRNATranscription',
                    'mRNADegradation',
                    'proteinTranslation',
                    'proteinDegradation',
                    'heavisideSigma',
                    'signalingTimescale',
                    'hillCoefficient',
                    'interactionStrength']
    boolodeDefaults = {
        'mRNATranscription':20.,
        'mRNADegradation':10.,
        'proteinTranslation':10.,
        'proteinDegradation':1.0,
        'heavisideSigma':10.,
        'signalingTimescale':5.0,
        'hillCoefficient':10.,
        'interactionStrength':1.0
    }
    kineticParameterDefaults = {}
    for rp in requiredPars:
        if rp not in parameterDefaults:
            print("%s is missing from the config file. Using default value.")            
        kineticParameterDefaults[rp] = parameterDefaults.get(rp, boolodeDefaults[rp])

    # Max level checks
    x_max = kineticParameterDefaults['mRNATranscription']/kineticParameterDefaults['mRNADegradation']
    y_max = x_max*(kineticParameterDefaults['proteinTranslation']/kineticParameterDefaults['proteinDegradation'])
    
    if x_max != 2.0:
        print('Warning: max(mRNA) != 2.0')
    if y_max != 20.:
        print('Warning: max(protein) != 20.0')

    return kineticParameterDefaults

def getRegulatorsInRule(rule, species, inputs):
    """
    Helper function to tokenize a rule into regulators.
    Returns three lists of regulator names.
    1. allreg is the list of all valid regulators
    2. regulatorySpecies are other model variables that are regulators
    3. inputreg are the regulators that are model inputs
    """
    rhs = rule
    rhs = rhs.replace('(',' ')
    rhs = rhs.replace(')',' ')
    tokens = rhs.split(' ')

    allreg = set([t for t in tokens if (t in species or t in inputs)])
    regulatorySpecies = set([t for t in tokens if t in species])
    inputreg = set([t for t in tokens if t in inputs])

    return((allreg, regulatorySpecies, inputreg))


def getSaneNval(size,lo=1.,hi=10.,mu=2.,sig=2.,identicalPars=False):
    """
    Generates a gaussian random number which is
    bounded by `lo` and `hi`

    :param size: number of random numbers to generate
    :type size: int
    :param lo: lower bound of sample range
    :type lo: float
    :param hi: upper bound of sample range
    :type hi: float
    :param mu: Mean of Gaussian distribution to be sampled from
    :type mu: float
    :param sigma: Standard deviation of Gausssian distribution to be sampled from
    :type sigma: float
    :param identicalPars: Flag to sample single value and return a list of identical values
    :type identicalPars: bool
    :returns:
        - K: list of sampled values
    """
    
    if identicalPars:
        k = np.random.normal(mu, sig)
        while k < lo or k > hi:
            k = np.random.normal(mu, sig)
        K = [k for i in range(size)]
    else:
        K = []
        for _ in range(size):
            k = np.random.normal(mu, sig)
            while k < lo or k > hi:
                k = np.random.normal(mu, sig)
            K.append(k)
    return K

def minmaxnorm(X):
    """Scales the values in X

    :param X: Input list of values to be scaled
    :type X: list
    :returns: 
        - N : list of values scaled between min and max values in input list
    """
    mix = min(X)
    mx = max(X)
    N = [(x-mix)/(mx-mix) for x in X]
    return N


def normalizeData(P):
    """
    Calls minmaxnorm() for each time series

    :returns:
        Pnorm : list of scaled values
    """
    Pnorm = []
    for i in range(np.shape(P)[1]):
        Pnorm.append(minmaxnorm(P[:,i]))
    return Pnorm


def normalizeExp(DF):
    """
    Calls minmaxnorm() for each gene across all experiments

    :param DF: Dataframe containing gene expression values as rows, time points as columns
    :type DF: pd.DataFrame
    :returns:
        newDF : DataFrame with values of DF scaled between min and max of DF
    """
    genes = DF.index
    newDF = DF.copy()
    Pnorm = []
    for g in genes:
        P = DF.loc[g].values
        newDF.loc[g] = minmaxnorm(P)
        
    return newDF

def get_ss(P):
    """
    Return the final time point of the simulated time course

    :param P: 2-D array containing time course of the model
    :type P: np.array
    :returns:
        - ss : np.array containing the last entry in each time series list 
    """
    ss = []
    ss = [p for p in P[-1,:]]
    return(ss)

def generateInputFiles(resultDF, BoolDF, withoutRules,
                       parameterInputsDF,tmax,numcells,
                       outPrefix=''):
    """
    Generates input files required from the Beeline pipeline

    :param resultDF: The simulation output, rows are genes, columns are "cells" or timepoints
    :type resultDF: pandas DataFrame
    :param outputfilenames: List of filenames generated containing individual time courses
    :type outputfilenames: list
    :param BoolDF: Dataframe containing rules
    :type BoolDF: pandas DataFrame
    :param withoutrules: List of nodes in input file without rules
    :type withoutrules: list
    :param parameterInputsPath: Specifies if there are any inputs to the model
    :type parameterInputsPath: str
    :param outPrefix: Prefix specifying target directory
    :type outPrefix: str (Optional)
    """
    
    print('1. refNetwork')
    refnet = []
    genes = set(BoolDF['Gene'].values)
    genes = genes.difference(set(withoutRules))
    inputs = withoutRules

    for g in genes:
        row = BoolDF[BoolDF['Gene'] == g]
        rhs = list(row['Rule'].values)[0]
        rule = list(row['Rule'].values)[0]
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')
        if len(withoutRules) == 0:
            inputs = []
            avoidthese = ['and','or', 'not', '']
        else:
            avoidthese = list(withoutRules)
            avoidthese.extend(['and','or', 'not', ''])

        regulators = [t for t in tokens if (t in genes or t in inputs) if t not in avoidthese]
        if 'not' in tokens:
            whereisnot = tokens.index('not')
        else:
            whereisnot = None
        for r in regulators:
            if whereisnot is None:
                ty = '+'
            else:
                if type(whereisnot) is int:
                    whereisnot = [whereisnot]
                
                if tokens.index(r) < whereisnot[0]:
                    ty = '+'
                else:
                    ty = '-'
             # Regulator is Gene1 and Target is Gene2
            refnet.append({'Gene2':g, 
                           'Gene1':r,
                           'Type':ty})
    refNetDF = pd.DataFrame(refnet)
    refNetDF.drop_duplicates(inplace=True)
    refNetDF.to_csv(str(outPrefix) + '/refNetwork.csv',sep=',',index=False)
    
    # PseudoTime.csv
    print('2. PseudoTime.csv')
    cellID = list(resultDF.columns)
    time = [float(c.split('_')[1].replace('-','.')) for c in cellID]
    experiment = [int(c.split('_')[0].split('E')[1]) for c in cellID]
    pseudotime = minmaxnorm(time)
    cellID = [c.replace('-','_') for c in cellID]

    PseudoTimeDict = {'Cell ID':cellID, 'PseudoTime':pseudotime,
                      'Time':time,'Experiment':experiment}
    PseudoTimeDF = pd.DataFrame(PseudoTimeDict)
    PseudoTimeDF.to_csv(str(outPrefix) + '/PseudoTime.csv',sep=',',index=False)
    PseudoTimeDF.index = PseudoTimeDF['Cell ID']
    
    # ExpressionData.csv
    if len(resultDF.columns) < 1e3:
        print('3. ExpressionData.csv')
        columns = list(resultDF.columns)
        columns = [c.replace('-','_') for c in columns]
        resultDF.columns = columns
        if parameterInputsDF is not None:
            resultDF = resultDF.drop(withoutRules, axis=0)
        resultDF.to_csv(str(outPrefix) + '/ExpressionData.csv',sep=',')
    else:
        print("Dataset too large."
              "\nSampling %d cells, one from each simulated trajectory." % numcells)
        times = np.random.choice([i for i in range(1,tmax*100)],numcells)
        expdf = pd.DataFrame(columns=['E' + str(i) + '_' + str(times[i])\
                                      for i in range(numcells)],
                             index=resultDF.index)
        for c in expdf.columns:
            expdf[c] = resultDF[c]
        expdf.to_csv(str(outPrefix) + '/ExpressionData.csv',sep=',')

def sampleTimeSeries(num_timepoints, expnum,\
                     tspan,  P,\
                     varmapper,timeIndex,
                     genelist, proteinlist,
                     header,writeProtein=False):
    """
    Returns pandas DataFrame with columns corresponding to 
    time points and rows corresponding to genes
    """
    revvarmapper = {v:k for k,v in varmapper.items()}
    experimentTimePoints = [h for h in header if 'E' + str(expnum) in h]
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    sampleDict = {}
    
    if writeProtein:
        # Write protein and mRNA to file
        for ri in varmapper.keys():
            sampleDict[varmapper[ri]] = {h:P[ri][ti] for h,ti in zip(experimentTimePoints,timeIndex)}
    else:
        # Default, only mRNA
        if len(proteinlist) == 0:
            for ri in rnaIndex:
                sampleDict[varmapper[ri]] = {h:P[ri][ti]\
                                             for h,ti in zip(experimentTimePoints,timeIndex)}
        else:
            speciesoi = [revvarmapper['p_' + p] for p in proteinlist]
            speciesoi.extend([revvarmapper['x_' + g] for g in genelist])
            result = pd.DataFrame(index=pd.Index([varmapper[i] for i in speciesoi]))

            for si in speciesoi:
                sampleDict[varmapper[si]] = {h:P[si][ti]\
                                             for h,ti in zip(experimentTimePoints,timeIndex)}

    sampleDF = pd.DataFrame(sampleDict)
    return(sampleDF)

def sampleCellFromTraj(cellid,
                       tspan,
                       P,
                       varmapper,sampleAt,
                       genelist, proteinlist,
                       header,writeProtein=False):
    """
    Returns pandas DataFrame with columns corresponding to 
    time points and rows corresponding to genes
    """
    revvarmapper = {v:k for k,v in varmapper.items()}
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    sampleDict = {}
    timepoint = int(header[cellid].split('_')[1])
    if writeProtein:
        # Write protein and mRNA to file
        for ri in varmapper.keys():
            sampleDict[varmapper[ri]] = {header[cellid]: P[ri][timepoint]}
    else:
        # Default, only mRNA
        if len(proteinlist) == 0:
            for ri in rnaIndex:
                sampleDict[varmapper[ri]] = {header[cellid]: P[ri][timepoint]}
        else:
            speciesoi = [revvarmapper['p_' + p] for p in proteinlist]
            speciesoi.extend([revvarmapper['x_' + g] for g in genelist])
            result = pd.DataFrame(index=pd.Index([varmapper[i] for i\
                                                  in speciesoi]))
            for si in speciesoi:
                sampleDict[varmapper[si]] = {header[cellid]: P[ri][timepoint]}

    sampleDF = pd.DataFrame(sampleDict)
    return(sampleDF)


def checkValidInputPath(path):
    """
    Returns dataframe of file at path.
    If path is not valid, returns empty dataframe.
    """
    df = pd.DataFrame()
    if path.is_file():
        df = pd.read_csv(path, sep='\t',engine='python')
    return(df)


def checkValidModelDefinitionPath(path, name):
    # Check if model defintion file exists
    if not os.path.isfile(path):
        print("Error in definition of job "  + name)
        print("Please specify path to Boolean model")
        return False
    return True
    
