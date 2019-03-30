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
#sys.path.insert(0, "/home/adyprat/anaconda3/envs/pyDSTool/lib/python3.6/site-packages/")
import PyDSTool as dst
from tqdm import tqdm

def readBooleanRules(path):
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
def createModelObject(ModelSpec,isStochastic):
    if isStochastic:
        tdomain = [0,200]
        times = dst.linspace(tdomain[0], tdomain[1],2000)
        xData = {'noise':np.random.randn(len(times))}
        my_input = dst.InterpolateTable({'tdata': times,
                                         'ics': xData,
                                         'name': 'interp'}).compute('interp')
        # plt.plot(my_input.sample()['noise'])
        # plt.show()
        # sys.exit()
        ModelDef = dst.args(name='bool',
                            varspecs=ModelSpec['varspecs'],
                            pars=ModelSpec['pars'],
                            ics=ModelSpec['ics'],
                            tdata=tdomain,
                            inputs={'noise':my_input.variables['noise']})
        DS = dst.Vode_ODEsystem(ModelDef)
    else:
        ModelDef = dst.args(name='bool',
                            varspecs=ModelSpec['varspecs'],
                            pars=ModelSpec['pars'],
                            ics=ModelSpec['ics'],
                            tdata=tdomain)
        DS = dst.Vode_ODEsystem(ModelDef)
        return DS

def simulateModel(DS):
    traj = DS.compute('test')
    points = traj.sample()
    return points

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
    par.update({'alpha_'+g:getSaneNval(mu=0.25,sig=1.,lo=0,hi=1) for g in genes})
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
                # TODO: For non random alphas to represent bool funcs modify this:
                par['a_' + currGen +'_'  + '_'.join(list(c))] = getSaneNval(mu=0.25,sig=1.,lo=0,hi=1) 
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
    with ('../' + outname,'w') as out:
        for k, v in ModelSpec['parameters']:
            out.write(k+'\t'+str(v))

def noise(x,t):

    return (1*np.sqrt(x))



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

        print(np.multiply(G(yn, tn),dWn))
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
    
def main():
    path = 'data/variables.txt'
    #path = 'data/test_vars.txt'
    DF = readBooleanRules(path)
    isStochastic = False
    ModelSpec = generateModelDict(DF)
    varmapper, parmapper = writeModelToFile(ModelSpec)
    import model
    variables = list(ModelSpec['ics'].keys())
    y0 = [ModelSpec['ics'][varmapper[i]] for i in range(len(varmapper.keys()))]
    
    rnaIndex = [i for i in range(len(varmapper.keys())) if 'x_' in varmapper[i]]
    
    proteinIndex = [i for i in range(len(variables)) if 'y_' in variables[i]]    
    pars = ModelSpec['pars'].values()
    t = np.linspace(0,100,1000)
    if not isStochastic:
        P = odeint(model.Model,y0,t,args=(pars,))
    else:
        P = eulersde(model.Model,noise,y0,t,pars)
        
    Pnorm = []
    # Min Max normalize
    for i in range(np.shape(P)[1]):
        Pnorm.append(minmaxnorm(P[:,i]))

        
    for ind in rnaIndex:
        plt.plot(t,Pnorm[ind],label=varmapper[ind])
        
    plt.legend()
    plt.show()
                        
if __name__ == "__main__":
    main()

# c in stoch sim is 0.05
