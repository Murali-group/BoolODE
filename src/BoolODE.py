#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'

import sys
import numpy as np
import scipy as sc
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
    withRules = list(DF['Gene'].values)
    allnodes = set()
    for ind,row in DF.iterrows():
        rhs = row['Rule']
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')

        reg = [t for t in tokens if t not in ['not','and','or','']]
        allnodes.update(set(reg))

    withoutRules = allnodes.difference(set(withRules))
    for n in withoutRules:
        DF = DF.append({'Gene':n,'Rule':n},ignore_index=True)
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

def generateModelDict(DF,identicalPars,samplePars):
    mRNATranscription = 20.
    mRNADegradation = mRNATranscription/2.
    proteinTranslation = 10
    proteinDegradation = 1
    hillThreshold = (proteinTranslation/proteinDegradation)*(mRNATranscription/mRNADegradation)/2
    hillCoefficient = 10
    genes = set(DF['Gene'].values)
    # Variables:
    varspecs = {'x_' + g:'' for g in genes}
    par = dict()
    if not samplePars:
        # Hill coefficients
        par.update({'n_'+g:hillCoefficient for g in genes})
        # Thresholds
        par.update({'k_'+g:hillThreshold for g in genes})
        # mRNA transcription rates
        par.update({'m_' + g:mRNATranscription for g in genes})
        # mRNA degradation rates
        par.update({'l_x_' + g:mRNADegradation for g in genes})
        # Protein translation rate
        par.update({'r_' + g:proteinTranslation for g in genes})
        # protein degradation rates
        par.update({'l_p_' + g:proteinDegradation for g in genes})
    else:
        # Hill coefficients
        Nvals = getSaneNval(len(genes),\
                        lo=0.5*hillCoefficient,\
                        hi=1.5*hillCoefficient,\
                        mu=hillCoefficient,\
                        sig=0.5*hillCoefficient,\
                        identicalPars=identicalPars)
        par.update({'n_'+g:n for g,n in zip(genes,Nvals)})
        # Thresholds
        kvals = getSaneNval(len(genes),\
                            lo=0.5*hillThreshold,
                            hi=1.5*hillThreshold,
                            mu=hillThreshold,
                            sig=0.5*hillThreshold,
                            identicalPars=identicalPars)
        par.update({'k_'+g:k for g,k in zip(genes,kvals)})
        # mRNA transcription rates
        mvals = getSaneNval(len(genes),\
                            lo=0.5*mRNATranscription,
                            hi=1.5*mRNATranscription,
                            mu=mRNATranscription,
                            sig=0.5*mRNATranscription,
                            identicalPars=identicalPars)
        par.update({'m_' + g:m for g,m in zip(genes,mvals)})
        # mRNA degradation rates
        lxvals = getSaneNval(len(genes),\
                             lo=0.5*mRNADegradation,
                             hi=1.5*mRNADegradation,
                             mu=mRNADegradation,
                             sig=0.5*mRNADegradation,
                             identicalPars=identicalPars)
        par.update({'l_x_' + g:lx for g,lx in zip(genes,lxvals)})
        # Protein translation rate
        rvals = getSaneNval(len(genes),\
                            lo=0.5*proteinTranslation,
                            hi=1.5*proteinTranslation,
                            mu=proteinTranslation,
                            sig=0.5*proteinTranslation,
                            identicalPars=identicalPars)
        par.update({'r_' + g:r for g,r in zip(genes,rvals)})
        # protein degradation rates
        lpvals = getSaneNval(len(genes),\
                             lo=0.5*proteinDegradation,
                             hi=1.5*proteinDegradation,
                             mu=proteinDegradation,
                             sig=0.5*proteinDegradation,
                             identicalPars=identicalPars)
        par.update({'l_p_' + g:lp for g,lp in zip(genes,lpvals)})
    
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
    xvals = [1. for _ in range(len(genes))]
    ics = {'x_' + g:x for g,x in zip(genes,xvals)}
    # This is reset in Experiment()
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
    c = 10.#4.
    return (c*np.sqrt(abs(x)))

def deltaW(N, m, h):
    # From sdeint implementation
    """Generate sequence of Wiener increments for m independent Wiener
    processes W_j(t) j=0..m-1 for each of N time intervals of length h.    
    Returns:
      dW (array of shape (N, m)): The [n, j] element has the value
      W_j((n+1)*h) - W_j(n*h) 
    """
    return np.random.normal(0.0, h, (N, m))

def eulersde(f,G,y0,tspan,pars,dW=None):
    # From sdeint implementation
    N = len(tspan)
    h = (tspan[N-1] - tspan[0])/(N - 1)
    maxtime = tspan[-1]
    # allocate space for result
    d = len(y0)
    y = np.zeros((N+1, d), dtype=type(y0[0]))

    if dW is None:
        # pre-generate Wiener increments (for d independent Wiener processes):
        dW = deltaW(N, d, h)
    y[0] = y0
    currtime = 0
    n = 0
   
    while currtime < maxtime:
        tn = currtime
        yn = y[n]
        dWn = dW[n,:]
        y[n+1] = yn + f(yn, tn,pars)*h + np.multiply(G(yn, tn),dWn)
        # Ensure positive terms
        for i in range(len(y[n+1])):
            if y[n+1][i] < 0:
                y[n+1][i] = yn[i]
        currtime += h
        n += 1 
    return y

def minmaxnorm(X):
    mix = min(X)
    mx = max(X)
    N = [(x-mix)/(mx-mix) for x in X]
    return N

def parseArgs(args):
    parser = OptionParser()
    parser.add_option('', '--max-time', type='int',default=20,
                      help='Total time of simulation. (Default = 20)')
    parser.add_option('', '--num-timepoints', type='int',default=100,
                      help='Number of time points to sample. (Default = 100)')
    parser.add_option('', '--num-experiments', type='int',default=30,
                      help='Number of experiments to perform. (Default = 30)')
    parser.add_option('-n', '--normalize-trajectory', action="store_true",default=False,
                      help="Min-Max normalize genes across all experiments")
    parser.add_option('-i', '--identical-pars', action="store_true",default=False,
                      help="Set single value to similar parameters\n"
                      "NOTE: Consider setting this if you are using --sample-pars.")
    parser.add_option('-s', '--sample-pars', action="store_true",default=False,
                      help="Sample rate parameters around the hardcoded means\n"
                      ", using 10% stand. dev.")
    parser.add_option('', '--write-protein', action="store_true",default=False,
                      help="Write both protein and RNA values to file. Useful for debugging.")
    parser.add_option('-b', '--burn-in', action="store_true",default=False,
                      help="Treats the first 25% of the time course as burnin\n"
                      ", samples from the rest.")        
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

def getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex, varmapper,revvarmapper):
    # Initialize
    new_ics = [0 for _ in range(len(varmapper.keys()))]
    # Set the mRNA ics
    for ind in rnaIndex:
        if ss[ind] < 0:
            ss[ind] = 0.0
        new_ics[ind] =  ss[ind]
        if new_ics[ind] < 0:
            new_ics[ind] = 0
    # Calculate the Protein ics based on mRNA levels
    proteinss = {}
    for i in rnaIndex:
        genename = varmapper[i].replace('x_','')
        proteinname = 'p_' + genename
        proteinss[proteinname] = ((ModelSpec['pars']['r_' + genename])/\
                                  (ModelSpec['pars']['l_p_' + genename]))\
                                  *new_ics[revvarmapper['x_' + genename]]
    # Finally record the protein ics values
    for ind in proteinIndex:
        new_ics[ind] = proteinss[varmapper[ind]]
    return(new_ics)

def Experiment(Model, ModelSpec,tspan, num_experiments,
               num_timepoints,
               varmapper, parmapper, outPrefix,
               burnin=False,writeProtein=False,
               normalizeTrajectory=False):
    pars = list(ModelSpec['pars'].values())
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    revvarmapper = {v:k for k,v in varmapper.items()}
    proteinIndex = [i for i in range(len(varmapper.keys())) if 'p_' in varmapper[i]]


    y0 = [ModelSpec['ics'][varmapper[i]] for i in range(len(varmapper.keys()))]
    
    ss = np.zeros(len(varmapper.keys()))
    for i,k in varmapper.items():
        if 'x_' in k:
            ss[i] = 1.

    outputfilenames = []
    for isStochastic in [True]: 
        # "WT" simulation
        result = pd.DataFrame(index=pd.Index([varmapper[i] for i in rnaIndex]))
        frames = []
        every = len(tspan)/num_timepoints
        if burnin:
            burninFraction = 0.25 # Chosen arbitrarily as 25%
            # Ignore the first 25% of the simulation
            startat = int(np.ceil(burninFraction*len(tspan)))
        else:
            startat = 0
            
        timeIndex = [i for i in range(startat,len(tspan)) if float(i)%float(every) == 0]
        header = ['E' + str(expnum) + '_' + str(round(tspan[tpoint],3)).replace('.','-')\
                  for expnum in range(num_experiments)\
                  for tpoint in timeIndex]
        
        for expnum in tqdm(range(num_experiments)):
            y0_exp = getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex, varmapper,revvarmapper)
            P = simulateModel(Model, y0_exp, pars, isStochastic, tspan)
            P = P.T
            # Extract Time points
            sampleDF = sampleTimeSeries(num_timepoints,expnum,\
                                        tspan, P,\
                                        varmapper,timeIndex, header,
                                        writeProtein=writeProtein)
            sampleDF = sampleDF.T
            frames.append(sampleDF)
        every = len(tspan)/num_timepoints

        result = pd.concat(frames,axis=1)
        result = result[header]
        indices = result.index
        newindices = [i.replace('x_','') for i in indices]
        result.index = pd.Index(newindices)

        if normalizeTrajectory:
            resultN = normalizeExp(result)
        else:
            resultN = result
            
        if isStochastic:
            name = 'stoch'
        else:
            name = 'ode'
        resultN.to_csv(outPrefix + name +'_experiment.txt',sep='\t')
        outputfilenames.append(outPrefix + name +'_experiment.txt')
    return outputfilenames
            
def sampleTimeSeries(num_timepoints, expnum,\
                     tspan,  P,\
                     varmapper,timeIndex,header,writeProtein=False):
    """
    Returns dictionary of dictionaries
    {gene : { exp_id : simulated value } }
    """
    every = int(len(tspan)/num_timepoints)
    experimentTimePoints = [h for h in header if 'E' + str(expnum) in h]
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    sampleDict = {}
    
    if writeProtein:
        # Write protein and mRNA to file
        for ri in varmapper.keys():
            sampleDict[varmapper[ri]] = {h:P[ri][ti] for h,ti in zip(experimentTimePoints,timeIndex)}
    else:
        # Default, only mRNA
        for ri in rnaIndex:
            sampleDict[varmapper[ri]] = {h:P[ri][ti] for h,ti in zip(experimentTimePoints,timeIndex)}

    sampleDF = pd.DataFrame(sampleDict)
    return(sampleDF)

def generateInputFiles(outputfilenames, BoolDF, outPrefix=''):
    for f in outputfilenames:
        syntheticDF = pd.read_csv(f,sep='\t',index_col=0)
        
        # ExpressionData.csv
        
        ExpDF = syntheticDF.copy()
        columns = list(ExpDF.columns)
        columns = [c.replace('-','_') for c in columns]
        ExpDF.columns = columns
        ExpDF.to_csv(outPrefix+'ExpressionData.csv',sep=',')
        
        # PseudoTime.csv
        cellID = list(syntheticDF.columns)
        time = [float(c.split('_')[1].replace('-','.')) for c in cellID]
        experiment = [int(c.split('_')[0].split('E')[1]) for c in cellID]
        pseudotime = minmaxnorm(time)
        cellID = [c.replace('-','_') for c in cellID]

        PseudoTimeDict = {'Cell ID':cellID, 'PseudoTime':pseudotime,
                          'Time':time,'Experiment':experiment}
        PseudoTimeDF = pd.DataFrame(PseudoTimeDict)
        PseudoTimeDF.to_csv(outPrefix + 'PseudoTime.csv',sep=',',index=False)

        # refnetwork
        refnet = []
        genes = set(BoolDF['Gene'].values)

        for g in genes:
            row = BoolDF[BoolDF['Gene'] == g]
            rhs = list(row['Rule'].values)[0]
            rule = list(row['Rule'].values)[0]
            rhs = rhs.replace('(',' ')
            rhs = rhs.replace(')',' ')
            tokens = rhs.split(' ')
            regulators = [t for t in tokens if t in genes if t not in ['and','or', 'not', '']]
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
                refnet.append({'Gene1':g,
                               'Gene2':r,
                               'Type':ty})
        refNetDF = pd.DataFrame(refnet)
        refNetDF.to_csv(outPrefix + 'refNetwork.csv',sep=',',index=False)
            
        

def main(args):
    opts, args = parseArgs(args)
    path = opts.path
    if path is None or len(path) == 0:
        print("Please specify path to Boolean model")
        sys.exit()

    tmax = opts.max_time
    numExperiments = opts.num_experiments
    numTimepoints = opts.num_timepoints
    identicalPars = opts.identical_pars
    burnin = opts.burn_in
    samplePars = opts.sample_pars
    outPrefix = opts.outPrefix

    writeProtein = opts.write_protein
    normalizeTrajectory = opts.normalize_trajectory
    
    timesteps = 100
    tspan = np.linspace(0,tmax,tmax*timesteps)
    
    DF = readBooleanRules(path)
    it = 0
    someexception = True
    while someexception:
        try:
            genesDict = {}
            
            ModelSpec = generateModelDict(DF,identicalPars,samplePars)
            varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
            parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
            writeModelToFile(ModelSpec)
            import model
            outputfilenames = Experiment(model.Model,ModelSpec,tspan,numExperiments,
                                         numTimepoints, varmapper, parmapper,
                                         outPrefix,burnin=burnin,
                                         writeProtein=writeProtein,
                                         normalizeTrajectory=normalizeTrajectory)
            generateInputFiles(outputfilenames,DF, outPrefix=outPrefix)
            print('Success!')
            
            someexception= False
            
        except FloatingPointError as e:
            it +=1 
            print(e,"\nattempt %d" %it)
        
    writeParametersToFile(ModelSpec)

if __name__ == "__main__":
    main(sys.argv)

