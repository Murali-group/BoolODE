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

def generateModelDict(DF,isStochastic):
    genes = set(DF['Gene'].values)
    # Variables:
    varspecs = {'x_' + g:'' for g in genes}

    par = {'n_'+g:getSaneNval() for g in genes}
    # Common Thresholds, uniform dist
    par.update({'k_'+g:getSaneNval(lo=0.01,hi=1.0,mu=0,sig=0.02) for g in genes})
    
    # basal activations, gaussian dist
    par.update({'a_'+g:getSaneNval(mu=0.25,sig=1.,lo=0,hi=1) for g in genes})
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
        Production = 'm_'+ currGen + '*(' +num + '/' + den + ')'
        Degradation = 'l_x_'  + currGen + '*x_' + currGen
        varspecs['x_' + currGen] =  Production \
                                   + '-' + Degradation \
                                   + stochasticTerm(isStochastic,Production,Degradation)\
                                   

    varspecs.update({'p_' + g:'r_'+g+'*'+'x_' +g + '- l_p_'+g+'*'+'p_' + g\
                     + stochasticTerm(isStochastic, 'r_'+g+'*'+'x_' +g,
                                      'l_p_'+g+'*'+'p_' + g)
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
        for i,v in enumerate(ModelSpec['varspecs'].keys()):
            out.write('    ' + v + ' = Y[' + str(i) + ']\n')
        
        for v,vdef in ModelSpec['varspecs'].iteritems():
            d =vdef.replace('noise','np.random.normal()')
            d = d.replace('^','**')
            out.write('    ' + 'd' + v + ' = ' + d  + '\n')
            outstr += 'd' + v + ', '
        for i,v in enumerate(ModelSpec['varspecs'].keys()):
            out.write('    if ' + v + ' < 0.0 :\n')
            # Reset to non zero value
            out.write('        ' + v + ' = Y['+str(i)+']\n')            

        out.write('    dY = [' + outstr+ ']\n')
        out.write('    return(dY)\n')
        out.write('#####################################################')
        
        
def stochasticTerm(isStochastic,Production, Degradation):
    c = 1 # 0.05From GNW 
    if isStochastic:
        # return(' + ' + str(c) \
        #        + '*noise*(('+Production+')^0.5 + ('+Degradation+')^0.5)')
        return(' + ' + str(c) \
               + '*noise*((abs('+Production+ ')+abs('+Degradation+'))^0.5)')        
    else:
        return('')
    

def writeParametersToFile(ModelSpec,outname='parameters.txt'):
    with ('../' + outname,'w') as out:
        for k, v in ModelSpec['parameters']:
            out.write(k+'\t'+str(v))


def rk4(yVector,currentTime,stepSize,Function,P):
    k1=np.array(Function(yVector,currentTime,P))
    V2=2.0*np.ones(len(k1))
    V3=6.0*np.ones(len(k1))
    k2=np.array(Function(yVector+0.5*np.ones(len(k1))*k1*stepSize,currentTime+stepSize*0.5,P))
    k3=np.array(Function(yVector+0.5*np.ones(len(k1))*k2*stepSize,currentTime+stepSize*0.5,P))
    k4=np.array(Function(yVector+k3*stepSize,currentTime+stepSize,P))
    yVectorNext=yVector+(k1+2.0*k2+2.0*k3+k4)*stepSize/6.0
    return yVectorNext

def solver(InitialCondition,TimeRange,stepSize,Function,P):
    startTime,stopTime=TimeRange
    Time=np.arange(startTime+stepSize,stopTime,stepSize)
    yVals=[InitialCondition]
    # for t in tqdm(Time):
    #     yVals.append(list(rk4(yVals[-1],t,stepSize,Function)))
    Y=InitialCondition
    t=startTime+stepSize
    PTIME=time.clock()
    for t in tqdm(Time):
        Y=list(rk4(yVals[-1],t,stepSize,Function,P))
        yVals.append(Y)
        #t=t+stepSize
        TOTALPTIME=time.clock()-PTIME
    return(yVals)

def main():
    path = 'data/variables.txt'
    # path = 'data/test_vars.txt'
    DF = readBooleanRules(path)
    isStochastic = True
    ModelSpec = generateModelDict(DF,isStochastic)
    writeModelToFile(ModelSpec)
    # t= np.linspace(0,100,1000)
    import model
    y0 = ModelSpec['ics'].values()
    pars = ModelSpec['pars'].values()
    #P = odeint(model.Model,y0,t,args=(pars,))
    
    tmin = 0
    tmax = 200.0
    stepsize = 0.01
    
    P,T = solver(y0, [tmin,tmax],stepsize,model.Model,pars)

    for i in range(0,len(P[0])):
        plt.plot(t,[p[i] for p in P])

    plt.show()
                        
if __name__ == "__main__":
    main()

# c in stoch sim is 0.05
