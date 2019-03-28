#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'

import sys
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations

# Uncomment on Aditya's machine
#sys.path.insert(0, "/home/adyprat/anaconda3/envs/pyDSTool/lib/python3.6/site-packages/")
import PyDSTool as dst

def readBooleanRules(path):
    DF = pd.read_csv("variables.txt",sep='\t')
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

def generateModelDict(DF):
    geneList = set(DF['Gene'].values)
    # Variables:
    varspecs = {'x_' + g:'' for g in genes}
    varspecs.update({'p_' + g:'r_'+g+'*'+'x_' +g + '- l_p_'+g+'*'+'p_' + g
                     for g in genes})

    par = {'n_'+g:getSaneNval() for g in genes}
    # Common Thresholds, uniform dist
    par.update({'k_'+g:getSaneNval(lo=0.01,hi=1.0,mu=0,sig=0.02) for g in genes})
    
    # basal activations, gaussian dist
    par.update({'a_'+g:getSaneNval(mu=0.25,sig=1.,lo=0,hi=1) for g in genes})
    par.update({'r_' + g:getSaneNval(lo=0.,hi=0.25,mu=0.0,sig=0.05) for g in genes})
    
    # mRNA degradation rates, currently identical for all mRNAs?
    #par.update({'l_x_' + g:np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100) for g in genes})
    lx = np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100)
    par.update({'l_x_' + g:lx for g in genes})

    # protein degradation rates, currently identical for all proteins?
    #par.update({'l_p_' + g:np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100) for g in genes})
    lp = np.log(2)/getSaneNval(mu=5,sig=10.,lo=0,hi=100)
    par.update({'l_p_' + g:lp for g in genes})

    
    par.update({'m_' + g:1.0 for g in genes})
    
    for i,row in DF.iterrows():
        rhs = row['Rule']
        rhs = rhs.replace('(',' ')
        rhs = rhs.replace(')',' ')
        tokens = rhs.split(' ')

        # This is currently unused
        # keywd = ['and','not', 'or', '']    
        reg = [t for t in tokens if t in genes]
        currGen = row['Gene']
        num = '( a_' + currGen
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
        varspecs['x_' + currGen] = 'l_x_'+ currGen + '*(' + num + '/' + den + '-' + 'x_' + currGen +')'
        
    # Initialize variables between 0 and 1, Doesn't matter.
    ics = {'x_' + g:getSaneNval(lo=0.,hi=1,mu=0.0,sig=1) for g in genes}
    ics.update({'p_' + g:ics['x_'+g]*par['r_'+g]/par['l_p_'+g] for g in genes})

    ModelSpec = {}
    ModelSpec['varspecs'] = varspecs
    ModelSpec['pars'] = pars
    ModelSpec['ics'] = ics
    # Default time
    ModelSpec['tdata'] = [0,200]
    return ModelSpec
    
def createModelObject(ModelSpec):
    ModelDef = dst.args(name='bool',
                        varspecs=ModelSpec['varspecs'],
                        pars=ModelSpec['par'],
                        ics=ModelSpec['ics'],
                        tdata=ModelSpec['tdata'])
    DS = dst.Vode_ODEsystem(ModelDef)
    return DS

def simulateModel(DS):
    traj = DS.compute('test')
    points = traj.sample()
    return points

def plotSimulation(points,varsToPlot):
    if len(varsToPlot) == 0:
        varsToPlot = list(points.keys())

    for v in varsToPlot:
        if v != 't':
            plt.plot(points['t'],points[v],label=v)
    plt.legend()
    plt.show()

def writeParametersToFile(ModelSpec,outname='parameters.txt'):
    with ('output/' + outname,'w') as out:
        for k, v in ModelSpec['parameters']:
            out.write(k+'\t'+str(v))

def main():
    path = '../data/variables.txt'
    DF = readBooleanRules(path)
    ModelSpec = generateModelDict(DF)
    DS = createModelObject(ModelSpec)
    points = simulateModel(DS)
    plotSimulation(points)
                        
if __name__ == "__main__":
    main()

