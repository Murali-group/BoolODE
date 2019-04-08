#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'

import sys
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations
from scipy.integrate import odeint
import time
# Uncomment on Aditya's machine
sys.path.insert(0, "/home/adyprat/anaconda3/envs/pyDSTool/lib/python3.6/site-packages/")
#import PyDSTool as dst
from tqdm import tqdm
from optparse import OptionParser

def readBooleanRules(path):
    print(path)
    DF = pd.read_csv(path,sep='\t')
    return DF

def getSaneNval(lo=1.,hi=10.,mu=2.,sig=2.):
    """
    Generates a gaussian random number which is
    bounded by lo and hi
    """
    k = np.random.normal(mu, sig)
    while k < lo or k > hi:
        k = np.random.normal(mu, sig)
    return k

def plotSimulation(points,varsToPlot=[]):
    if len(varsToPlot) == 0:
        varsToPlot = list(points.keys())

    for v in varsToPlot:
        if v != 't':
            plt.plot(points['t'],points[v],label=v)
    plt.legend()
    plt.show()

def generateModelDict(DF):
    genes = set(DF['Gene'].values)
    # Variables:
    varspecs = {'x_' + g:'' for g in genes}

    par = {'n_'+g:getSaneNval() for g in genes}
    # Common Thresholds, uniform dist
    par.update({'k_'+g:getSaneNval(lo=0.01,hi=1.0,mu=0,sig=0.02) for g in genes})
    
    # basal activations, gaussian dist
    par.update({'r_' + g:getSaneNval(lo=0.,hi=0.25,mu=0.0,sig=0.05) for g in genes})
    
    # mRNA degradation rates, currently identical for all mRNAs?
    par.update({'l_x_' + g:np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100) for g in genes})
    #lx = np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100)
    #par.update({'l_x_' + g:lx for g in genes})

    # protein degradation rates, currently identical for all proteins?
    par.update({'l_p_' + g:np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100) for g in genes})
    #lp = np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100)
    #par.update({'l_p_' + g:lp for g in genes})

    par.update({'m_' + g:getSaneNval(mu=5,sig=10.,lo=0,hi=100) for g in genes})
    

    genespace = {}
    for i,row in DF.iterrows():
        tempStr = row['Gene'] + " = 0"  
        exec(tempStr, genespace)

    
    
    for i,row in DF.iterrows():
        # Basal alpha
        #exec('booleval = 1', genespace) 
        #print(genespace['booleval'])
        #print(booleval)
        exec('booleval = ' + row['Rule'], genespace) 
        par['alpha_'+row['Gene']] = int(genespace['booleval'])#getSaneNval(mu=0.25,sig=1.,lo=0,hi=1) for g in genes})
    
    for i,row in DF.iterrows():
        rhs = row['Rule']
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')

        # This is currently unused
        # keywd = ['and','not', 'or', '']    
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
                #print(genespace['a'],genespace['b'],genespace['c']) 

                for geneInList in c:
                    exec(geneInList + ' = 1', genespace)
                #print(genespace['a'],genespace['b'],genespace['c']) 
                #print(row['Rule'])
                exec('boolval = ' + row['Rule'], genespace)

                par['a_' + currGen +'_'  + '_'.join(list(c))] = int(genespace['boolval']) #getSaneNval(mu=0.25,sig=1.,lo=0,hi=1)
                #print(par['a_'+currGen+'_'+'_'.join(list(c))], genespace['boolval'])
        num += ' )'
        den += ' )'
        #varspecs['x_' + currGen] = 'l_x_'+ currGen + '*' + num + '/' + den + '-' + 'l_x_'  + currGen +'*x_' + currGen
        Production = 'm_'+ currGen + '*(' +num + '/' + den + ')'
        Degradation = 'l_x_'  + currGen + '*x_' + currGen
        varspecs['x_' + currGen] =  Production \
                                   + '-' + Degradation #\
                                   # + stochasticTerm(isStochastic,Production,Degradation)\
                                   

    varspecs.update({'p_' + g:'r_'+g+'*'+'x_' +g + '- l_p_'+g+'*'+'p_' + g\
                     # + stochasticTerm(isStochastic, 'r_'+g+'*'+'x_' +g,
                     #                  'l_p_'+g+'*'+'p_' + g)
                     for g in genes})
        
    # Initialize variables between 0 and 1, Doesn't matter.
    ics = {'x_' + g:getSaneNval(lo=0.,hi=1,mu=0.0,sig=1) for g in genes}
    ics.update({'p_' + g:ics['x_'+g]*par['r_'+g]/par['l_p_'+g] for g in genes})

    ModelSpec = {}
    ModelSpec['varspecs'] = varspecs
    ModelSpec['pars'] = par
    ModelSpec['ics'] = ics
    # Default time
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
    c = 0.05
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
    parser.add_option('', '--max-time', type='int',
                      help='Total time of simulation')
    parser.add_option('', '--num-timepoints', type='int',
                      help='Number of time points to sample')
    parser.add_option('', '--num-experiments', type='int',
                      help='Number of experiments to perform')
    parser.add_option('', '--outPrefix', type = 'str',
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
            ss[ind] = 1
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
               varmapper, parmapper, outPrefix = ''):
    pars = list(ModelSpec['pars'].values())
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    revvarmapper = {v:k for k,v in varmapper.items()}
    proteinIndex = [i for i in range(len(varmapper.keys())) if 'p_' in varmapper[i]]
    
    y0 = [ModelSpec['ics'][varmapper[i]] for i in range(len(varmapper.keys()))]
    # First do the ODE simulations, no noise, then stoch
    for isStochastic in [True,False]: 
        # "WT" simulation\
        #result = None
        result = pd.DataFrame(index=pd.Index([varmapper[i] for i in rnaIndex]))
        frames = []
        
        ss = get_ss(simulateModel(Model, y0, pars, isStochastic, tspan))
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
    
def plotRNA(Pnorm, rnaIndex):
    # Visualize
    for ind in rnaIndex:
        plt.plot(t,Pnorm[ind],label=varmapper[ind])
        
    plt.legend()
    plt.show()
    
def main(args):
    opts, args = parseArgs(args)
    
    path = opts.path
    tmax = opts.max_time
    num_experiments = opts.num_experiments
    num_timepoints = opts.num_timepoints
    DF = readBooleanRules(path)
    
    genesDict = {}
       	

    ModelSpec = generateModelDict(DF)
    
    writeParametersToFile(ModelSpec)
    
    varmapper, parmapper = writeModelToFile(ModelSpec)
    
    import model
    
    tspan = np.linspace(0,tmax,tmax*10)

    Experiment(model.Model,ModelSpec,tspan,num_experiments,num_timepoints, varmapper,parmapper, opts.outPrefix) 
                        
if __name__ == "__main__":
    main(sys.argv)

# c in stoch sim is 0.05
