import numpy as np
import pandas as pd
import os
import sys

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

def writeModelToFile(ModelSpec, prefix=''):
    """
    Writes model to file as a python function, which is then imported.

    :param ModelSpec: ODE equations stored in a dictionary
    :type ModelSpec: dict
    :param prefix: Optional argument that specifies prefix to model filename
    :type prefix: str
    :returns:
        - dir_path : path to output directory
    :rtype: str
        
    """
    varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
    parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}
    dir_path = os.path.dirname(os.path.realpath(__file__))    
    print("I am in " + dir_path)
    with open(dir_path+"/"+prefix+'model.py','w') as out:
        out.write('#####################################################\n')
        out.write('import numpy as np\n')
        out.write('# This file is created automatically\n')
        out.write('def Model(Y,t,pars):\n')
        out.write('    # Parameters\n')
        par_names = sorted(ModelSpec['pars'].keys())
        for i,p in enumerate(par_names):
            out.write('    ' + p + ' = pars[' + str(i) + ']\n')
        outstr = ''
        out.write('    # Variables\n')
        for i in range(len(varmapper.keys())):
            out.write('    ' + varmapper[i] + ' = Y[' + str(i) + ']\n')
            outstr += 'd' + varmapper[i] + ','
        for i in range(len(varmapper.keys())):
            vdef = ModelSpec['varspecs'][varmapper[i]]
            vdef = vdef.replace('^','**')
            out.write('    d' + varmapper[i] + ' = '+vdef+'\n')
            
        out.write('    dY = np.array([' + outstr+ '])\n')
        out.write('    return(dY)\n')
        out.write('#####################################################')
    return dir_path 

def writeParametersToFile(ModelSpec, outPrefix, outname='parameters.txt'):
    """
    Writes dictionary of parameters to file

    :param ModelSpec: Model definition dictionary containing the keys 'ics', 'pars', and 'varspec' 
    :type ModelSpec: dict 
    :param outPrefix: Prefix to output folder name
    :type outPrefix: str
    :param outname: User defined name for parameters file. Default is parameters.txt
    :type outname: str (Optional)
    """
    with open(str(outPrefix) + '/'  + outname,'w') as out:
        for k, v in ModelSpec['pars'].items():
            out.write(k+'\t'+str(v) + '\n')

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

def generateInputFiles(resultDF, outputfilenames, BoolDF, withoutRules,
                       parameterInputsDF,
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
    
    for f in outputfilenames:
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
        if len(resultDF.columns) < 1e5:
            print('3. ExpressionData.csv')
            columns = list(resultDF.columns)
            columns = [c.replace('-','_') for c in columns]
            resultDF.columns = columns
            if parameterInputsDF is not None:
                resultDF = resultDF.drop(withoutRules, axis=0)
            resultDF.to_csv(str(outPrefix) + '/ExpressionData.csv',sep=',')
        else:
            print('Dataset too large. Skipping generation of ExpressionData.csv.\n Please sample from simulations.')
            

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
    #experimentTimePoints = [h for h in header if 'E' + str(expnum) in h]
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
    if os.path.isfile(path):
        df = pd.read_csv(path,sep='\t')
        return(df)
    else:
        return(None)

def checkValidModelDefinitionPath(path, name):
    # Check if model defintion file exists
    if not os.path.isfile(path):
        print("Error in definition of job "  + name)
        print("Please specify path to Boolean model")
        return False
    return True
    
