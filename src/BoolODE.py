#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'

import sys
import numpy as np
import scipy as sc
#import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations
from scipy.integrate import odeint
import time
import importlib
import warnings
# Uncomment on Aditya's machine
sys.path.insert(0, "/home/adyprat/anaconda3/envs/pyDSTool/lib/python3.6/site-packages/")
from tqdm import tqdm
from optparse import OptionParser

np.seterr(all='raise')

def readBooleanRules(path):
    DF = pd.read_csv(path,sep='\t')
    return DF

def getSaneNval(size,lo=1.,hi=10.,mu=2.,sig=2.,identicalPars=False):
    """
    Generates a gaussian random number which is
    bounded by lo and hi
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

def generateModelDict(DF,identicalPars):
    genes = set(DF['Gene'].values)
    # Variables:
    varspecs = {'x_' + g:'' for g in genes}
    # Hill coefficients
    N = getSaneNval(len(genes),identicalPars=identicalPars)
    par = {'n_'+g:n for g,n in zip(genes,N)}
    # Thresholds
    kvals = getSaneNval(len(genes),identicalPars=identicalPars,lo=0.01,hi=1.0,mu=0.5,sig=0.02)
    par.update({'k_'+g:k for g,k in zip(genes,kvals)})
    # mRNA transcription rates
    mvals = getSaneNval(len(genes),identicalPars=identicalPars,mu=10,sig=5.,lo=0,hi=100)
    par.update({'m_' + g:m for g,m in zip(genes,mvals)})
    # mRNA degradation rates
    lxvals = getSaneNval(len(genes),identicalPars=identicalPars,mu=10,sig=5.,lo=0,hi=100)
    par.update({'l_x_' + g:np.log(2)/lx for g,lx in zip(genes,lxvals)})
    # Protein translation rate
    rvals = getSaneNval(len(genes),identicalPars=identicalPars,lo=50.,hi=100.0,mu=75.0,sig=10)
    par.update({'r_' + g:r for g,r in zip(genes,rvals)})
    # protein degradation rates
    lpvals = getSaneNval(len(genes),identicalPars=identicalPars,mu=10,sig=5.,lo=0,hi=100)
    par.update({'l_p_' + g:np.log(2)/lp for g,lp in zip(genes,lpvals)})
    
    # Initialize new namespace
    genespace = {}
    for i,row in DF.iterrows():
        # Initialize variables to 0
        tempStr = row['Gene'] + " = 0"  
        exec(tempStr, genespace)
    
    for i,row in DF.iterrows():
        # Basal alpha:
        # Execute the rule to figure out
        # the value of alpha
        exec('booleval = ' + row['Rule'], genespace) 
        par['alpha_'+row['Gene']] = int(genespace['booleval'])
    
    for i,row in DF.iterrows():
        rhs = row['Rule']
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')

        reg = [t for t in tokens if t in genes]
        currGen = row['Gene']
        num = '( alpha_' + currGen
        den = '( 1'
        for i in range(1,len(reg) + 1):
            for c in combinations(reg,i):
                # Create the hill function terms for each regulator
                hills = ['(p_'+ci+'/k_'+ci+')^n_'+ci for ci in c]
                mult = '*'.join(hills)
                # Create Numerator and Denominator
                den += ' +' +  mult                
                num += ' + a_' + currGen +'_'  + '_'.join(list(c)) + '*' + mult
                
                for i1, row1 in DF.iterrows():
                    exec(row1['Gene'] + ' = 0', genespace)

                for geneInList in c:
                    exec(geneInList + ' = 1', genespace)
                    
                exec('boolval = ' + row['Rule'], genespace)

                par['a_' + currGen +'_'  + '_'.join(list(c))] = int(genespace['boolval']) 

        num += ' )'
        den += ' )'
        
        Production = 'm_'+ currGen + '*(' +num + '/' + den + ')'
        Degradation = 'l_x_'  + currGen + '*x_' + currGen
        varspecs['x_' + currGen] =  Production \
                                   + '-' + Degradation
                                   
    varspecs.update({'p_' + g:'r_'+g+'*'+'x_' +g + '- l_p_'+g+'*'+'p_' + g\
                     for g in genes})
        
    # Initialize variables between 0 and 1, Doesn't matter.
    xvals = getSaneNval(len(genes),identicalPars=identicalPars,lo=0.,hi=100,mu=10.0,sig=5)
    ics = {'x_' + g:x for g,x in zip(genes,xvals)}
    ics.update({'p_' + g:ics['x_'+g]*par['r_'+g]/par['l_p_'+g] for g in genes})

    ModelSpec = {}
    ModelSpec['varspecs'] = varspecs
    ModelSpec['pars'] = par
    ModelSpec['ics'] = ics
    return ModelSpec

def writeModelToFile(ModelSpec):
    varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
    parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
    
    with open('src/model.py','w') as out:
        out.write('#####################################################\n')
        out.write('import numpy as np\n')
        out.write('# This file is created automatically\n')
        out.write('def Model(Y,t,pars):\n')
        out.write('    # Parameters\n')
        for i,p in enumerate(ModelSpec['pars'].keys()):
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
    return varmapper,parmapper
        
def writeParametersToFile(ModelSpec,outname='parameters.txt'):
    with open(outname,'w') as out:
        for k, v in ModelSpec['pars'].items():
            out.write(k+'\t'+str(v) + '\n')

def noise(x,t):
    c = 1e-2
    return (c*np.sqrt(x))

def deltaW(N, m, h):
    # From sdeint implementation
    """Generate sequence of Wiener increments for m independent Wiener
    processes W_j(t) j=0..m-1 for each of N time intervals of length h.    
    Returns:
      dW (array of shape (N, m)): The [n, j] element has the value
      W_j((n+1)*h) - W_j(n*h) 
    """
    return np.random.normal(0.0, np.sqrt(h), (N, m))

def eulersde(f,G,y0,tspan,pars,dW=None):
    # From sdeint implementation
    N = len(tspan)
    h = (tspan[N-1] - tspan[0])/(N - 1)
    # allocate space for result
    d = len(y0)
    y = np.zeros((N, d), dtype=type(y0[0]))

    if dW is None:
        # pre-generate Wiener increments (for d independent Wiener processes):
        dW = deltaW(N - 1, d, h)
    y[0] = y0
    for n in range(0, N-1):
        tn = tspan[n]
        yn = y[n]
        dWn = dW[n,:]

        y[n+1] = yn + f(yn, tn,pars)*h + np.multiply(G(yn, tn),dWn)
        for i in range(len(y[n+1])):
            if y[n+1][i] < 0:
                y[n+1][i] = yn[i]
    return y

def minmaxnorm(X):
    mix = min(X)
    mx = max(X)
    N = [(x-mix)/(mx-mix) for x in X]
    return N

def parseArgs(args):
    parser = OptionParser()
    parser.add_option('', '--max-time', type='int',default=100,
                      help='Total time of simulation')
    parser.add_option('', '--num-timepoints', type='int',default=10,
                      help='Number of time points to sample')
    parser.add_option('', '--num-experiments', type='int',default=10,
                      help='Number of experiments to perform')
    parser.add_option('-i', '--identical-pars', action="store_true",default=False,
                      help='Set single value to similar parameters')    
    parser.add_option('', '--outPrefix', type = 'str',default='',
                      help='Prefix for output files.')
    parser.add_option('', '--path', type='str',
                      help='Path to boolean model file')    
    (opts, args) = parser.parse_args(args)

    return opts, args

def simulateModel(Model, y0, parameters,isStochastic, tspan):
    if not isStochastic:
        P = odeint(Model,y0,tspan,args=(parameters,))
    else:
        P = eulersde(Model,noise,y0,tspan,parameters)
    return(P)

def normalizeData(P):
    """
    Calls minmaxnorm() for each time series
    """
    Pnorm = []
    for i in range(np.shape(P)[1]):
        Pnorm.append(minmaxnorm(P[:,i]))
    return Pnorm

def get_ss(P):
    """
    Returns the last entry in each time series list 
    """
    ss = []
    ss = [p for p in P[-1,:]]
    return(ss)

def getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex, varmapper,revvarmapper):
    # Initialize
    new_ics = [0 for _ in range(len(varmapper.keys()))]
    # Set the mRNA ics
    for ind in rnaIndex:
        if ss[ind] < 0:
            ss[ind] = 1#1
        new_ics[ind] =  np.random.normal(ss[ind],ss[ind]*0.25)
        if new_ics[ind] < 0:
            new_ics[ind] = 0
    # Calculate the Protein ics based on mRNA levels
    proteinss = {}
    for i in rnaIndex:
        genename = varmapper[i].replace('x_','')
        proteinname = 'p_' + genename
        proteinss[proteinname] = ((ModelSpec['pars']['r_' + genename])/(ModelSpec['pars']['l_p_' + genename]))*new_ics[revvarmapper['x_' + genename]]
    # Finally record the protein ics values
    for ind in proteinIndex:
        new_ics[ind] = proteinss[varmapper[ind]]
    return(new_ics)

def Experiment(Model, ModelSpec,tspan, num_experiments,
               num_timepoints,
               varmapper, parmapper, outPrefix):
    pars = list(ModelSpec['pars'].values())
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    revvarmapper = {v:k for k,v in varmapper.items()}
    proteinIndex = [i for i in range(len(varmapper.keys())) if 'p_' in varmapper[i]]
    
    y0 = [ModelSpec['ics'][varmapper[i]] for i in range(len(varmapper.keys()))]
    
    # First do an ODE simulation to get ss
    ss = get_ss(simulateModel(Model, y0, pars, False, np.linspace(0,200,1000)))

    for isStochastic in [True,False]: 
        # "WT" simulation
        result = pd.DataFrame(index=pd.Index([varmapper[i] for i in rnaIndex]))
        frames = []
        
        for expnum in range(num_experiments):
            y0_exp = getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex, varmapper,revvarmapper)
            P = simulateModel(Model, y0_exp, pars, isStochastic, tspan)
            # Min Max normalize; Do this last!
            Pnorm = normalizeData(P)
            # Extract Time points
            sampleDF = sampleTimeSeries(num_timepoints,expnum,tspan,rnaIndex,Pnorm,varmapper)
            sampleDF = sampleDF.T
            frames.append(sampleDF)
        every = len(tspan)/num_timepoints
        timeIndex = [i for i in range(1,len(tspan)) if i%every == 0]
        columns = []
        for expnum in range(num_experiments):
            for tpoint in timeIndex:
                columns.append('E' + str(expnum) + '_' + str(int(tspan[tpoint])))
        result = pd.concat(frames,axis=1)
        result = result[columns]
        if isStochastic:
            name = 'stoch'
        else:
            name = 'ode'
        result.to_csv(outPrefix + name +'_experiment.txt',sep='\t')
            
def sampleTimeSeries(num_timepoints,expnum,tspan,rnaIndex,P, varmapper):
    every = len(tspan)/num_timepoints
    timeIndex = [i for i in range(1,len(tspan)) if i%every == 0]
    sampleDict = {}
    for ri in rnaIndex:
        sampleDict[varmapper[ri]] = {'E'+str(expnum)+'_'+str(int(tspan[ti])):\
                                     P[ri][ti] for ti in timeIndex}
    sampleDF = pd.DataFrame(sampleDict)
    return(sampleDF)
    
def main(args):
    opts, args = parseArgs(args)
    path = opts.path
    if path is None or len(path) == 0:
        print("Please specify path to Boolean model")
        sys.exit()

    tmax = opts.max_time
    num_experiments = opts.num_experiments
    num_timepoints = opts.num_timepoints
    tspan = np.linspace(0,tmax,tmax*10)
    identicalPars = opts.identical_pars
    DF = readBooleanRules(path)
    it = 0
    someexception = True
    while someexception:
        try:
            genesDict = {}
            
            ModelSpec = generateModelDict(DF,identicalPars)
            varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
            parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
            writeModelToFile(ModelSpec)
            import model
            Experiment(model.Model,ModelSpec,tspan,num_experiments,
                       num_timepoints, varmapper, parmapper,
                       opts.outPrefix)
            print('Success!')
            someexception= False
        except FloatingPointError as e:
            it +=1 
            print(e,"\nattempt %d" %it)
        
    writeParametersToFile(ModelSpec)

if __name__ == "__main__":
    main(sys.argv)

