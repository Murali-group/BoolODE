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
import ast
np.seterr(all='raise')
import os

def readBooleanRules(path, parameterInputsPath):
    """
    Parameters
    ----------
    path : str
        Path to Boolean Rule file
    parameterInputsPath : str
        Path to file containing input parameters to the model
    Returns
    -------
    DF : pandas DataFrame
        Dataframe containing rules
    withoutrules : list
        list of nodes in input file without rules
    """
    DF = pd.read_csv(path,sep='\t',engine='python')
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
        if len(parameterInputsPath) == 0:
            print(n, "has no rule, adding self-activation.")
            DF = DF.append({'Gene':n,'Rule':n},ignore_index=True)
        else:
            print("Treating %s as parameter" % {n})
    return DF, withoutRules

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

def getParameters(DF,identicalPars,
                  samplePars,
                  genelist,
                  proteinlist,
                  withoutRules,
                  parameterInputsDF,
                  interactionStrengthDF):
    """
    Create dictionary of parameters and values. Assigns
    parameter values by evaluating the Boolean expression
    for each variable.
    Parameters:
    -----------
    DF : pandas DataFrame
        Table of values with two columns, 'Gene' specifies
        target, 'Rule' specifies Boolean function
    identicalPars : bool
        Passed to getSaneNval to set identical parameters
    samplePars : bool
        Sample kinetic parameters using a Gaussian distribution 
        centered around the default parameters
    parameterInputsDF : pandas DataFrame
        Optional table that specifies parameter values. Useful to specify
        experimental conditions.
    interactionStrengthDF : pandas DataFrame
        Optional table that specifies interaction strengths. When
        not specified, default strength is set to 1.
    Returns:
    --------
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
    hillThreshold = (proteinTranslation/proteinDegradation)*\
                    (mRNATranscription/mRNADegradation)/2
    y_max = (proteinTranslation/proteinDegradation)*\
            (mRNATranscription/mRNADegradation)
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
    ## Whether to sample random parameters or stick to defaults
    if samplePars:
        print("Sampling parameters")
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            sampledParameterValues = getSaneNval(len(species),\
                                                 lo=0.5*parDefault,\
                                                 hi=1.5*parDefault,\
                                                 mu=parDefault,\
                                                 sig=0.5*parDefault,\
                                                 identicalPars=identicalPars)
            for sp, sparval in zip(species, sampledParameterValues):
                if sp in genelist:
                    par[parPrefix + sp] = sparval                
            
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            sampledParameterValues = getSaneNval(len(species),\
                                                 lo=0.5*parDefault,\
                                                 hi=1.5*parDefault,\
                                                 mu=parDefault,\
                                                 sig=0.5*parDefault,\
                                                 identicalPars=identicalPars)
            for sp, sparval in zip(species, sampledParameterValues):
                if sp in genelist:
                    par[parPrefix + sp] = sparval                
        
    else:
        print("Fixing rate parameters to defaults")
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            for sp in species:
                par[parPrefix + sp] = parDefault

        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            for sp in species:
                if sp in genelist:
                    par[parPrefix + sp] = parDefault

    return par, parameterInputsDict

def generateModelDict(DF,identicalPars,
                      samplePars,
                      withoutRules,
                      speciesTypeDF,
                      parameterInputsDF,
                      interactionStrengthDF):
    """
    Take a DataFrame object with Boolean rules,
    construct ODE equations for each variable.
    Takes optional parameter inputs and interaction
    strengths
    Parameters:
    ----------
    DF : pandas DataFrame
        Table of values with two columns, 'Gene' specifies
        target, 'Rule' specifies Boolean function
    identicalPars : bool
        Passed to getSaneNval to set identical parameters
    samplePars : bool
        Sample kinetic parameters using a Gaussian distribution 
        centered around the default parameters
    speciesTypeDF : pandas DataFrame
        Table defining type of each species. Takes values 
        'gene' or 'protein'
    parameterInputsDF : pandas DataFrame
        Optional table that specifies parameter values. Useful to specify
        experimental conditions.
    interactionStrengthDF : pandas DataFrame
        Optional table that specifies interaction strengths. When
        not specified, default strength is set to 1.
    Returns:
    --------
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

                
    print(genelist)
    print(proteinlist)

    par, parameterInputs = getParameters(DF,identicalPars,
                                         samplePars,
                                         genelist,
                                         proteinlist,
                                         withoutRules,
                                         parameterInputsDF,
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

    for i,row in tqdm(DF.iterrows()):
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
    # varspecs.update({'p_' + g:'r_'+g+'*'+'x_' +g + '- l_p_'+g+'*'+'p_' + g\
    #                  for g in species})
        
    # Initialize variables between 0 and 1, Doesn't matter.
    xvals = [1. for _ in range(len(genelist))]
    pvals = [20. for _ in range(len(proteinlist))]    
    ics = {}

    for sp,xv in zip(genelist,xvals):
        ics['x_' + sp] = xv
        ics['p_' + sp] = 0
    for sp, pv in zip(proteinlist,pvals):
        ics['p_' + sp] = pv

    ModelSpec = {}
    ModelSpec['varspecs'] = varspecs
    ModelSpec['pars'] = par
    ModelSpec['ics'] = ics
    return ModelSpec, parameterInputs, genelist, proteinlist

def writeModelToFile(ModelSpec, prefix=''):
    varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
    parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
    
    with open('src/'+prefix+'model.py','w') as out:
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
    # Controls noise proportional to
    # square root of activity
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
    parser.add_option('', '--inputs', type='str',default='',
                      help='Path to input parameter files')    
    parser.add_option('', '--ics', type='str',default='',
                      help='Path to list of initial conditions')
    parser.add_option('', '--strengths', type='str',default='',
                      help='Path to list of interaction strengths')
    parser.add_option('', '--species-type', type='str',default='',
                      help="Path to list of molecular species type file."\
                      "Useful to specify proteins/genes")
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

def getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex,
                        genelist, proteinlist,
                        varmapper,revvarmapper):

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

def Experiment(Model, ModelSpec,tspan, num_experiments,
               num_timepoints,
               varmapper, parmapper,
               genelist, proteinlist,
               outPrefix,icsDF,
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
            
            # if 'x_' in k:
            #     if k.replace('x_','') in genes:
            #         ss[i] = icsmap[k.replace('x_','')]
            #     else:
            #         ss[i] = 0.01
        

    outputfilenames = []
    for isStochastic in [True]: 
        # "WT" simulation
        if len(proteinlist) == 0:
            result = pd.DataFrame(index=pd.Index([varmapper[i] for i in rnaIndex]))
        else:
            speciesoi = [revvarmapper['p_' + p] for p in proteinlist]
            speciesoi.extend([revvarmapper['x_' + g] for g in genelist])
            result = pd.DataFrame(index=pd.Index([varmapper[i] for i in speciesoi]))
            
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
            y0_exp = getInitialCondition(ss, ModelSpec, rnaIndex, proteinIndex,
                                         genelist, proteinlist,
                                         varmapper,revvarmapper)

            P = simulateModel(Model, y0_exp, pars, isStochastic, tspan)
            P = P.T
            # Extract Time points
            sampleDF = sampleTimeSeries(num_timepoints,expnum,\
                                        tspan, P,\
                                        varmapper, timeIndex,
                                        genelist, proteinlist,                                        
                                        header,
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

        # if output directory doesn't exist
        # create one!
        if len(outPrefix) > 0:
            if '/' in outPrefix:
                outDir = '/'.join(outPrefix.split('/')[:-1])
                if not os.path.exists(outDir):
                    print(outDir, "does not exist, creating it...")
                    os.makedirs(outDir)

        resultN.to_csv(outPrefix + name +'_experiment.txt',sep='\t')
        outputfilenames.append(outPrefix + name +'_experiment.txt')
    return outputfilenames
            
def sampleTimeSeries(num_timepoints, expnum,\
                     tspan,  P,\
                     varmapper,timeIndex,
                     genelist, proteinlist,
                     header,writeProtein=False):
    """
    Returns dictionary of dictionaries
    {gene : { exp_id : simulated value } }
    """
    revvarmapper = {v:k for k,v in varmapper.items()}
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

def generateInputFiles(outputfilenames, BoolDF, withoutRules,
                       parameterInputsPath,
                       outPrefix=''):
    for f in outputfilenames:
        syntheticDF = pd.read_csv(f,sep='\t',index_col=0,engine='python')
        
        # ExpressionData.csv

        
        ExpDF = syntheticDF.copy()
        columns = list(ExpDF.columns)
        columns = [c.replace('-','_') for c in columns]
        ExpDF.columns = columns
        if len(parameterInputsPath) == 0:
            ExpDF = ExpDF.drop(withoutRules, axis=0)
        ExpDF.to_csv(outPrefix+'ExpressionData.csv',sep=',')
        ExpDF.drop([col for col in ExpDF.columns if ExpDF[col].max()\
                    < 0.1*ExpDF.values.max()],axis=1,inplace=True)
        ExpDF.to_csv(outPrefix+'ExpressionData-dropped.csv',sep=',')        
        
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
        PseudoTimeDF.index = PseudoTimeDF['Cell ID']
        PseudoTimeDF.loc[ExpDF.columns].to_csv(outPrefix +\
                                               'PseudoTime-dropped.csv', sep = ',', index = False)
        # refnetwork
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
    parameterInputsPath = opts.inputs
    icsPath = opts.ics    
    writeProtein = opts.write_protein
    normalizeTrajectory = opts.normalize_trajectory
    interactionStrengthPath = opts.strengths
    speciesTypePath = opts.species_type
    
    if len(parameterInputsPath) > 0: 
        parameterInputsDF = pd.read_csv(parameterInputsPath,sep='\t')
    else:
        parameterInputsDF = None

    if len(icsPath) > 0: 
        icsDF = pd.read_csv(icsPath,sep='\t',engine='python')
    else:
        icsDF = None

    if len(interactionStrengthPath) > 0:
        print("Interaction Strengths are supplied")
        interactionStrengthDF = pd.read_csv(interactionStrengthPath,sep='\t',engine='python')
    else:
        interactionStrengthDF = None

    if len(speciesTypePath) > 0: 
        speciesTypeDF = pd.read_csv(speciesTypePath,sep='\t')
    else:
        speciesTypeDF = None
    
    timesteps = 100
    tspan = np.linspace(0,tmax,tmax*timesteps)
    
    DF, withoutRules = readBooleanRules(path, parameterInputsPath)
    if len(withoutRules) == 0:
        withoutRules = []
    it = 0
    someexception = True
    while someexception:
        try:
            genesDict = {}
            
            ModelSpec,\
                parameterInputs,\
                genelist,\
                proteinlist = generateModelDict(DF,identicalPars,
                                                           samplePars,
                                                           withoutRules,
                                                           speciesTypeDF,
                                                           parameterInputsDF,
                                                           interactionStrengthDF)
            # FIXME : ask user to pass row if parameter input file
            # Hardcoded. We only care about one input vector, say the first one
            if len(parameterInputsPath) > 0:
                ModelSpec['pars'].update(parameterInputs[0])
            varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
            parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}    
            writeModelToFile(ModelSpec)
            import model
            outputfilenames = Experiment(model.Model,ModelSpec,tspan,numExperiments,
                                         numTimepoints, varmapper, parmapper,
                                         genelist,proteinlist,
                                         outPrefix, icsDF,
                                         burnin=burnin,
                                         writeProtein=writeProtein,
                                         normalizeTrajectory=normalizeTrajectory)
            generateInputFiles(outputfilenames,DF,
                               withoutRules,
                               parameterInputsPath,
                               outPrefix=outPrefix)
            print('Success!')
            
            someexception= False
            
        except FloatingPointError as e:
            it +=1 
            print(e,"\nattempt %d" %it)
        
    writeParametersToFile(ModelSpec)

if __name__ == "__main__":
    main(sys.argv)

