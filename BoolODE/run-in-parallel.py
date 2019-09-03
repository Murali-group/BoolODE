from multiprocessing import Process, Lock
import pandas as pd
import BoolODE as bo
import os
import importlib.machinery
import numpy as np

models = [
    [{'model':'dyn-linear.txt',
     'input':'',
      'ics':'dyn-linear_ics.txt',
      'strengths':'',#'dyn-linear_strengths.txt'
}],
    [{'model':'dyn-bifurcating.txt',
      'input':'',
      'ics':'dyn-bifurcating_ics.txt',
      'strengths':'',#'dyn-bifurcating_strengths.txt'
}],
    [{'model':'dyn-linear-long.txt',
      'input':'',
      'ics':'dyn-linear-long_ics.txt',
      'strengths':'',#'dyn-linear-long_strengths.txt'
}],
    [{'model':'dyn-cycle.txt',
      'input':'',
      'ics':'dyn-cycle_ics.txt',
      'strengths':'',#'dyn-cycle_strengths.txt'
}],
    [{'model':'dyn-consecutive-bifurcating.txt',
      'input':'',
      'ics':'dyn-consecutive-bifurcating_ics.txt',
      'strengths':'',#'dyn-consecutive-bifurcating_strengths.txt'
}],
    [{'model':'dyn-bifurcating-converging.txt',
      'input':'',
      'ics':'dyn-bifurcating-converging_ics.txt',
      'strengths':'',#'dyn-bifurcating-converging_strengths.txt'
}],
    [{'model':'dyn-trifurcating.txt',
      'input':'',
      'ics':'dyn-trifurcating_ics.txt',
      'strengths':'',#'dyn-trifurcating_strengths.txt'
}],
    [{'model':'dyn-bifurcating-loop.txt',
      'input':'',
      'ics':'dyn-bifurcating-loop_ics.txt',
      'strengths':'',#'dyn-bifurcating-loop_strengths.txt'
}]
]

def startSimulations(argdict):

    tmax = 15.0
    numExperiments = 100
    numTimepoints = 100
    identicalPars = False
    burnin= False
    writeProtein= False
    normalizeTrajectory= False
    samplePars = False
    timesteps = 100                               
    tspan = np.linspace(0,tmax,tmax*timesteps)
    prefix = 'data/'
    path = prefix + argdict['model']
    outPrefix = 'dyn-output/' + argdict['model'].split('.txt')[0] + '_'
    interactionStrengthPath = prefix + argdict['strengths']


    if len(argdict['strengths']) == 0:
        interactionStrengthPath = ''
    else:
        interactionStrengthPath = prefix + argdict['strengths']
    
    if len(interactionStrengthPath) > 0:
        interactionStrengthDF = pd.read_csv(interactionStrengthPath,sep='\t',engine='python')
    else:
        interactionStrengthDF = None
        
    if len(argdict['input']) == 0:
        parameterInputsPath = ''
    else:
        parameterInputsPath = prefix + argdict['input']
    if len(argdict['ics']) == 0:
        icsPath = ''
    else:
        icsPath = prefix + argdict['ics']

    if len(parameterInputsPath) > 0: 
        parameterInputsDF = pd.read_csv(parameterInputsPath,sep='\t',engine='python')
    else:
        parameterInputsDF = None
        
    if len(icsPath) > 0: 
        icsDF = pd.read_csv(icsPath,sep='\t',engine='python')
    else:
        icsDF = None

    DF, withoutRules = bo.readBooleanRules(path, parameterInputsPath)

    ModelSpec, parameterInputs = bo.generateModelDict(DF,identicalPars,
                                                      samplePars,
                                                      withoutRules,
                                                      parameterInputsDF,
                                                      interactionStrengthDF)
    
    # Hardcoded. We only care about one input vector, say the first one
    if len(parameterInputsPath) > 0:
        ModelSpec['pars'].update(parameterInputs[0])
    varmapper = {i:var for i,var in enumerate(ModelSpec['varspecs'].keys())}
    parmapper = {i:par for i,par in enumerate(ModelSpec['pars'].keys())}
    dirname = argdict['model'].replace('.txt','')
    if not os.path.exists('src/'+ dirname + '/'):
        os.mkdir('src/'+ dirname + '/')
    
    bo.writeModelToFile(ModelSpec,prefix=dirname + '/')

    model = importlib.machinery.SourceFileLoader('model', 'src/' +dirname+'/model.py').load_module()
    outputfilenames = bo.Experiment(model.Model,ModelSpec,tspan,numExperiments,
                                 numTimepoints, varmapper, parmapper,
                                 outPrefix, icsDF,
                                 burnin=burnin,
                                 writeProtein=writeProtein,
                                 normalizeTrajectory=normalizeTrajectory)
    
    bo.generateInputFiles(outputfilenames,DF,
                       withoutRules,
                       parameterInputsPath,
                       outPrefix=outPrefix)

if __name__ == '__main__':
    lock = Lock()
    for i in range(len(models)):
        print(models[i][0]['model'])
        Process(target=startSimulations,args=(models[i])).start()
    
