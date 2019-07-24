import numpy as np
import pandas as pd
import os
import sys

def getSaneNval(size,lo=1.,hi=10.,mu=2.,sig=2.,identicalPars=False):
    """
    Generates a gaussian random number which is
    bounded by lo and hi
    Parameters
    ----------
    size : int
        number of random numbers to generate
    lo : float
        lower bound of sample range
    hi : float
        upper bound of sample range
    mu : float
        Mean of Gaussian distribution to be 
        sampled from
    sigma : float
        Standard deviation of Gausssian distribution
        to be sampled from
    identicalPars : bool
        Flag to sample single value and return a 
        list of identical values
    Returns 
    -------
    K : list
        list of sampled values
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
    return dir_path # varmapper,parmapper, 

def writeParametersToFile(ModelSpec, outPrefix, outname='parameters.txt'):
    with open(outPrefix + outname,'w') as out:
        for k, v in ModelSpec['pars'].items():
            out.write(k+'\t'+str(v) + '\n')

def minmaxnorm(X):
    mix = min(X)
    mx = max(X)
    N = [(x-mix)/(mx-mix) for x in X]
    return N


def normalizeData(P):
    """
    Calls minmaxnorm() for each time series
    """
    Pnorm = []
    for i in range(np.shape(P)[1]):
        Pnorm.append(minmaxnorm(P[:,i]))
    return Pnorm


def normalizeExp(DF):
    """
    Calls minmaxnorm() for each gene across all experiments
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
    Returns the last entry in each time series list 
    """
    ss = []
    ss = [p for p in P[-1,:]]
    return(ss)

def generateInputFiles(syntheticDF, outputfilenames, BoolDF, withoutRules,
                       parameterInputsPath,
                       outPrefix=''):
    for f in outputfilenames:
        #syntheticDF = pd.read_csv(f,sep='\t',index_col=0,engine='python')
        
        # refnetwork
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
        refNetDF.to_csv(outPrefix + 'refNetwork.csv',sep=',',index=False)
        
        # PseudoTime.csv
        print('2. PseudoTime.csv')
        cellID = list(syntheticDF.columns)
        time = [float(c.split('_')[1].replace('-','.')) for c in cellID]
        experiment = [int(c.split('_')[0].split('E')[1]) for c in cellID]
        pseudotime = minmaxnorm(time)
        cellID = [c.replace('-','_') for c in cellID]

        PseudoTimeDict = {'Cell ID':cellID, 'PseudoTime':pseudotime,
                          'Time':time,'Experiment':experiment}
        PseudoTimeDF = pd.DataFrame(PseudoTimeDict)
        PseudoTimeDF.to_csv(outPrefix + 'PseudoTime.csv',sep=',',index=False)
        PseudoTimeDF.index = PseudoTimeDF['Cell ID']
        PseudoTimeDF.loc[cellID].to_csv(outPrefix +\
                                               'PseudoTime-dropped.csv', sep = ',', index = False)
        
        # ExpressionData.csv
        print('3. ExpressionData.csv')
        columns = list(syntheticDF.columns)
        columns = [c.replace('-','_') for c in columns]
        syntheticDF.columns = columns
        if len(parameterInputsPath) == 0:
            syntheticDF = syntheticDF.drop(withoutRules, axis=0)
        syntheticDF.to_csv(outPrefix+'ExpressionData.csv',sep=',')
        

