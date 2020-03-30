#!/usr/bin/env python
# coding: utf-8
__author__ = 'Amogh Jalihal'
import os
import sys
import ast
import yaml
import time
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
import multiprocessing as mp
from optparse import OptionParser
from itertools import combinations
from scipy.integrate import odeint
from sklearn.cluster import KMeans
from typing import Dict, List
# local imports
from BoolODE import utils
from BoolODE import simulator 
from importlib.machinery import SourceFileLoader

class GenerateModel:
    """
    Generate an ODE model from the Boolean rule file specified by the user.
    
    Parameters:
    -----------
    :param settings: User defined job settings. This is created in boolode.do_post_processing()
    :type settings: dict
    :param parameterInputsDF: User specified parmeter inputs. Specified in config YAML as a path.
    :type parameterInputsDF: pandas DataFrame
    :param parameterSetDF: User specified kinetic parameters. Specified in config YAML as a path.
    :type parameterSetDF: pandas DataFrame
    :param interactionStrengthDF: User specified interaction strengths. Specified in config YAML as a path
    :type interactionStrengthDF: pandas DataFrame
    :param withRules: Nodes/genes in Boolean rules file with an associated Boolean rule
    :type withRules: list
    :param withoutRules: Nodes/genes in Boolean rules file that do not have an associated Boolean rule. If found, such nodes will be treated as ??
    :type withoutRules: list
    """
    def __init__(self, settings, parameterInputsDF, parameterSetDF, interactionStrengthDF) -> None:
        """
        initialize.
        """
        self.settings = settings
        self.parameterInputsDF = parameterInputsDF
        self.parameterSetDF = parameterSetDF
        self.interactionStrengthDF = interactionStrengthDF
        # Class variables
        self.withRules = list()
        self.withoutRules = list()
        self.allnodes = set()
        self.varspecs = dict()
        self.varmapper = dict()
        self.par = dict()
        self.parmapper = dict()
        self.kineticParameterDefaults = dict()
        self.inputs = list()
        self.genelist = list()
        self.proteinlist = list()
        self.nodeTypeDF = pd.DataFrame()
        self.df = pd.DataFrame()
        self.ModelSpec = dict()
        self.path_to_ode_model = str()
        # Read the model definition
        # 1. populate self.df
        # 2. store genelist, withRules, withoutRules, allnodes
        # 3. initialize varspecs
        self.readBooleanRules()
        # Create model parameters and assign values
        self.getParameters()
        # Create the model dictionary
        self.generateModelDict()
        # Write ODE model to file
        self.writeModelToFile()
        # Write parameters to file
        self.writeParametersToFile()
        
    def readBooleanRules(self):
        """
        Reads a rule file from path and stores the rules in a dictionary
        
    
        :returns:
            - df: Dataframe containing rules
            - withoutrules: list of nodes in input file without rules
        """
        self.df = pd.read_csv(self.settings['modelpath'], sep='\t', engine='python')
        
        self.withRules = list(self.df['Gene'].values)

        ## Extract all nodes from the boolean rules
        ## Note that not all nodes might have a rule attached to them
        for ind,row in self.df.iterrows():
            rhs = row['Rule']
            rhs = rhs.replace('(',' ')
            rhs = rhs.replace(')',' ')
            tokens = rhs.split(' ')
            reg = [t for t in tokens if t not in ['not','and','or','']]
            self.allnodes.update(set(reg))
    
        self.withoutRules = list(self.allnodes.difference(set(self.withRules)))

        ## Every node without a rule is treated as follows:
        ## If the user has specified a Parameter Input file treat as parameter, else 
        for n in self.withoutRules:
            if not self.parameterInputsDF.empty\
               and n in self.parameterInputsDF['Input'] :
                print("Treating %s as parameter" % {n})                
            else:
                print(n, "has no rule, adding self-activation.")
                self.df = self.df.append({'Gene':n,'Rule':n}, ignore_index=True)
                self.withRules.append(n)
                self.withoutRules.remove(n)
                
        if self.settings['add_dummy']:
            self.addDummyGenes()
        
        # Variables:
        ## Check:
        ## If nodeType is provided, initialize the
        ## correct node names.
        ## Rule:
        ## Every gene (x_node) gets a corresponding
        ## protein (p_node). The gene equation contains
        ## regulatory logic, protein equation contains production
        ## degradation terms.
        ## Every 'protein' (p_node) only gets one variable,
        ## and the corresponding equation contains the regulatory
        ## logic term. There is no production/degradation associated
        ## with this node.

        # Create empty equations
        if self.nodeTypeDF.empty:
            # Assume everything is a gene, so make the corresponding protein
            for node in self.withRules:
                self.varspecs['x_' + node] = ''
                self.varspecs['p_' + node] = ''
                self.genelist.append(node)
        else:
            # TODO Consider removing support for node type.
            for i, row in self.nodeTypeDF.iterrows():
                node = row['Node']
                if row['Type'] == 'protein':
                    self.varspecs.update({'p_' + node:''})
                    self.proteinlist.append(node)
                elif row['Type'] == 'gene':
                    self.varspecs.update({'x_' + node: '',
                                          'p_' + node:''})
                    self.genelist.append(node)
                else:
                    print(row['Type'],"not recognized. Treating as gene.")
                    self.varspecs.update({'x_' + node: '',
                                          'p_' + node: ''})
                    self.genelist.append(node)
    
    def addDummyGenes(self):
        ## Add dummy genes
        ## determinisitcally add parents
        print('using max_parents=' + str(self.settings['max_parents']))
        num_dummy = 2*len(allnodes)
        for dg in range(num_dummy):
            ## Variable number of parents
            # DF = DF.append({'Gene':'dummy' + str(dg),
            #                 'Rule':' or '.join([s for s in \
            #                                     np.random.choice(list(allnodes),
            #                                                      size=np.random.choice([i+1 for i in range(len(allnodes)-1)]))])},
            #                ignore_index = True)
            self.df = self.df.append({'Gene':'dummy' + str(dg),
                                      'Rule':' or '.join([s for s in \
                                                          np.random.choice(list(allnodes),
                                                                           size=max_parents)])},
                                     ignore_index = True)
        self.df.to_csv(self.setting['outPrefix'] + 'rules-with-added-genes.csv')        
        
    def getParameters(self):
        """
        Create dictionary of parameters and values. Assigns
        parameter values by evaluating the Boolean expression
        for each variable.
    
            pars : dict
                Dictionary of parameters 
        """
    
        self.kineticParameterDefaults = utils.loadParameterValues()
        
        ## Derived parameters:
        ## The threshold is calculated as the half max value of the species
        x_max = self.kineticParameterDefaults['mRNATranscription']/self.kineticParameterDefaults['mRNADegradation']
        y_max = x_max*(self.kineticParameterDefaults['proteinTranslation']/self.kineticParameterDefaults['proteinDegradation'])
        
        hillThreshold = y_max/2
        heavisideOmega = 2./y_max

        self.kineticParameterDefaults['x_max'] = x_max
        self.kineticParameterDefaults['y_max'] = y_max
        self.kineticParameterDefaults['hillThreshold'] = hillThreshold
        self.kineticParameterDefaults['heavisideOmega'] = heavisideOmega        
    
        parameterNamePrefixAndDefaultsAll = {
            # Hill coefficients
            'n_':self.kineticParameterDefaults['hillCoefficient'],
            # Thresholds
            'k_':hillThreshold,
            # Soft-Heaviside function steepness
            'sigmaH_':self.kineticParameterDefaults['heavisideSigma']
        }     
        parameterNamePrefixAndDefaultsGenes = {
            # mRNA transcription rates
            'm_':self.kineticParameterDefaults['mRNATranscription'],
            # mRNA degradation rates
            'l_x_':self.kineticParameterDefaults['mRNADegradation'],
            # Protein translation rate
            'r_':self.kineticParameterDefaults['proteinTranslation'],
            # protein degradation rates
            'l_p_':self.kineticParameterDefaults['proteinDegradation']
        }     
        
        par = dict()
        ## If there is even one signaling protein,
        ## create the y_max parameter
        if len(self.proteinlist) > 0:
            self.par['y_max'] = y_max
            self.par['signalingtimescale'] = signalingtimescale
            
        interactionStrengths = {}
    
        ## Carry out a series of checks to handle user defined input to model

        ## Check 1: handle "input parameters"
        if not self.parameterInputsDF.empty:
            # TODO The name 'parameter inputs' is confusing/misleading, considering renaming?
            self.addParameterInputs()

        ## Check 2: include user-specified interaction strengths
        if not self.interactionStrengthDF.empty:
            self.addInteractionStrengths(hillThreshold)

        ## Assigining values to each kinetic parameter
        ## Create a parameter of each parameter type for each gene.
        ## If samplePars=True, sample values from normal distributions,
        ## else, set all parameters to defaults
        if self.settings['sample_pars']:
            self.assignSampledParameterValues(parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes)
        else:
            self.assignDefaultParameterValues(parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes)

        # Final check:
        # Update user defined parameter values
        if not self.parameterSetDF.empty:
            for pname, pval in self.parameterSetDF.iterrows():
                self.par[pname] = pvalue[1]
                
    def addParameterInputs(self):
        """
        Check for user specified input parameters.
        If a parameter input file is specified,
        Every row in parameterInputsDF contains two
        columns:
        - 'Input' specifies the name of input parameter
        - 'Value' specified the input state, typically 1 or 0, but
          for HIGH and LOW inputs respectively, but any float between
          0 and 1 is valid input.
        Now, we create a fake hill term for this input, even
        though there is no equation corresponding to it. Thus,
        we create a hillThreshold and hillCoefficient 
        term for it. This is useful so we can treat this input
        parameter as a "regulator" and leave the structure of
        the logic term unchanged.
        TODO Add appropriate parameters if heaviside
        """
        parameterInputsDict  = {}
        hillThreshold = self.kineticParameterDefaults['hillThreshold']
        hillCoefficient = self.kineticParameterDefaults['hillCoefficient']
        for i, row in self.parameterInputsDF.iterrows():
            # initialize
            #parameterInputsDict[i] = {p:0 for p in self.withoutRules} ## Suspicious
            inputParam = row['Input']            
            if inputParam not in self.withoutRules:
                print(inputParam, "not found nodes without rules. Skipping")
            else:

                ## This can take values between 0 and 1
                if row['Value'] > 1.0:
                    inputValue = 1.0
                elif row['Value'] < 0.0:
                    inputValue = 0.0
                else:
                    inputValue = row['Value']
                # The following line sets the appropriate constant input
                # Note: 2*hillThreshold is the max value that a protein can take
                parameterInputsDict[inputParam] = inputValue*2*hillThreshold
        # Add input parameters, set to the value calculated above
        self.par.update({p:v for p, v in parameterInputsDict.items()})
        if self.settings['modeltype'] == 'hill':
            self.par.update({'k_'+p:hillThreshold for p in parameterInputsDict.keys()})
            self.par.update({'n_'+p:hillCoefficient for p in parameterInputsDict.keys()})

    def addInteractionStrengths(self, hillThreshold):
        ## Check 2:
        ## If interaction strengths are specified.
        ## Create a new parameter specifying the strength
        ## of a given interaction. For all other equations in
        ## which "regulator" appears in, its hillThreshold value
        ## will remain the default value.
        ## TODO Add appropriate parameters if heaviside                    
        for i,row in self.interactionStrengthDF.iterrows():
            regulator = row['Gene2']
            target = row['Gene1']
            interactionStrength = row['Strength']
            if self.settings['modeltype'] == 'hill':
                self.par.update({'k_' + regulator + '_' + target: hillThreshold/interactionStrength})
            elif self.settings['modeltype'] == 'heaviside':
                self.par.update({'w_' + regulator + '_' + target: hillThreshold/interactionStrength})
        

    def assignDefaultParameterValues(self,
                                     parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes):
        print("Fixing rate parameters to defaults")
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            for node in self.withRules:
                self.par[parPrefix + node] = parDefault
                
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            for node in self.withRules:
                if node in self.genelist:
                    self.par[parPrefix + node] = parDefault

    def assignSampledParameterValues(self,
                                     parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes):
        """ 
        Sample kinetic parameters from a normal distribution, but in a range +- 10% of 
        the default value.
        """
        print("Sampling parameter values")
        print("Using std=" + str(self.settings['sample_std']))
        lomult = 0.9
        himult = 1.1
        # Sample hillthreshold and hillcoeff
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            sampledParameterValues = utils.getSaneNval(len(self.withRules),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=self.settings['sample_std']*parDefault,\
                                                 identicalPars=self.settings['identical_pars'])
            for node, sparval in zip(self.withRules, sampledParameterValues):
                if node in self.genelist:
                    self.par[parPrefix + node] = sparval
                    
        transcriptionRate = 0.0 # Explain
        mRNADegradationRate = 0.0 # Explain
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            sampledParameterValues = utils.getSaneNval(len(self.withRules),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=self.settings['sample_std']*parDefault,\
                                                 identicalPars=self.settings['identical_pars'])
            # # This is very unclear! Why are we doing this?
            # if identicalPars:
            #     if parPrefix == 'm_':
            #         transcriptionRate = sampledParameterValues[0]
            #     if parPrefix == 'l_x_':
            #         mRNADegradationRate = sampledParameterValues[0]
            # else:
            #     if parPrefix == 'm_':
            #         transcriptionRate = parDefault
            #     if parPrefix == 'l_x_':
            #         mRNADegradationRate = parDefault
            
            for node, sparval in zip(self.withRules, sampledParameterValues):
                if node in self.genelist:
                    self.par[parPrefix + node] = sparval

    def createRegulatoryTerms(self, currgene, combinationOfRegulators,
                              strengthSpecified,
                              regulatorsWithStrength,
                              regSpecies):
        if self.settings['modeltype'] == 'hill':
            # Create the hill function terms for each regulator
            hills = []
            for reg in combinationOfRegulators:
                if strengthSpecified and (reg in regulatorsWithStrength):
                    hillThresholdName = 'k_' + reg+ '_' + currgene
                else:
                    hillThresholdName = 'k_' + reg

                if reg in regSpecies:
                    hills.append('(p_'+ reg +'/'+hillThresholdName+')^n_'+ reg)
                else:
                    # Note: Only proteins can be regulators
                    hills.append('('+ reg +'/'+hillThresholdName+')^n_'+ reg)

            mult = '*'.join(hills)
            return mult
        elif self.settings['modeltype'] == 'heaviside':
            terms = []
            for reg in combinationOfRegulators:
                terms.append('p_' + reg)
            mult = '*'.join(terms)
            return mult        
                    
    def generateModelDict(self):
        """
        Take a DataFrame object with Boolean rules,
        construct ODE equations for each variable.
        Takes optional parameter inputs and interaction
        strengths
    
        :param DF: Table of values with two columns, 'Gene' specifies target, 'Rule' specifies Boolean function
        :type DF: pandas DataFrame
        :param identicalPars: Passed to utils.getSaneNval to set identical parameters
        :type identicalPars: bool
        :param samplePars: Sample kinetic parameters using a Gaussian distribution  centered around the default parameters
        :type samplePars: bool
        :param sampleStd: Sample from a distribution of mean=par, sigma=sampleStd*par
        :type sampleStd: float
        :param speciesTypeDF: Table defining type of each species. Takes values  'gene' or 'protein'
        :type speciesTypeDF: pandas DataFrame
        :param parameterInputsDF: Optional table that specifies parameter input values. Useful to specify experimental conditions.
        :type parameterInputsDF: pandas DataFrame
        :param parameterSetDF: Optional table that specifies a predefined parameter set.
        :type parameterSetDF: pandas DataFrame
        :param interactionStrengthDF: Optional table that specifies interaction strengths. When not specified, default strength is set to 1.
        :type interactionStrengthDF: pandas DataFrame
        :param modeltype: Specify type of modeling framework. ['hill', 'heaviside']
        :type modeltype: str
        :returns:
            ModelSpec : dict
                Dictionary of dictionaries. 
                'varspecs' {variable:ODE equation}
                'ics' {variable:initial condition}
                'pars' {parameter:parameter value}
        """
                    
        ##########################################################
        ## Assign values to alpha parameters representing the logic
        ## relationships between variables
        # Initialize new namespace
        boolodespace = {}

        # Initialize empty equations
        for node in self.withRules:
            # Initialize species to 0
            tempStr = node + " = 0"  
            exec(tempStr, boolodespace)

        # If there is no rule correspondng to a node, it is either
        # a user specified parameter input, or it is assigned a self loop.
        if not self.parameterInputsDF.empty:
            self.inputs = set(self.withoutRules)
            for k in self.parameterInputsDF['Input']:
                # Initialize variables to 0
                tempStr = k + " = 0"  
                exec(tempStr, boolodespace)

        for i,row in self.df.iterrows():
            # Basal expression:
            # Execute the rule to figure out
            # the value of alpha_0 or omega_0
            exec('booleval = ' + row['Rule'], boolodespace)         
            if self.settings['modeltype'] == 'hill':
                self.par['alpha_'+row['Gene']] = int(boolodespace['booleval'])
            elif self.settings['modeltype'] == 'heaviside':
                self.par['omega_' + row['Gene']] = utils.heavisideThreshold(boolodespace['booleval'])
        # End initialization

        # Assign values to logic parameters, and construct expressions
        for i,row in self.df.iterrows():
            ## Parse Boolean rule to get list of regulators
            allreg, regSpecies, regInputs = utils.getRegulatorsInRule(row['Rule'],
                                                                      self.withRules,
                                                                      self.inputs)

            # Basal expression term
            currSp = row['Gene']
            if self.settings['modeltype'] == 'hill':
                num = '( alpha_' + currSp
                den = '( 1'
            elif self.settings['modeltype'] == 'heaviside':
               exponent = '- sigmaH_' + currSp +'*( omega_' + currSp
    
            strengthSpecified = False
            regulatorsWithStrength = set()
            if not self.interactionStrengthDF.empty:
                if currSp in self.interactionStrengthDF['Gene1'].values:
                    strengthSpecified = True
                    # Get the list of regulators of currSp
                    # whose interaction strengths have been specified
                    regulatorsWithStrength = set(self.interactionStrengthDF[\
                                                                       self.interactionStrengthDF['Gene1']\
                                                                       == currSp]['Gene2'].values)
                    print(regulatorsWithStrength)
                    
            # Loop over combinations of regulators        
            for i in range(1,len(allreg) + 1):
                for combinationOfRegulators in combinations(allreg,i):
                    regulatorExpression = self.createRegulatoryTerms(currSp, combinationOfRegulators,
                                                                      strengthSpecified,
                                                                      regulatorsWithStrength,
                                                                      regSpecies)
                    if self.settings['modeltype'] == 'hill':
                        # Create Numerator and Denominator
                        den += ' +' +  regulatorExpression
                        num += ' + a_' + currSp +'_'  + '_'.join(list(combinationOfRegulators)) + '*' + regulatorExpression
                    elif self.settings['modeltype'] == 'heaviside':
                        exponent += ' + w_' + currSp + '_' + '_'.join(list(combinationOfRegulators)) +'*' + regulatorExpression
    
                    # evaluate rule to assign values to parameters
                    ##################################################
                    for node in self.withRules:                 #
                        exec(node + ' = 0', boolodespace)  #                        
                    # Set each regulator to ON, evaluate rule. we are looping
                    # over all such combinations of regulators
                    for geneInList in combinationOfRegulators:     #
                        exec(geneInList + ' = 1', boolodespace)    #
                                                                   #
                    exec('boolval = ' + row['Rule'], boolodespace) #
                    ##################################################
    
                    if self.settings['modeltype'] == 'hill':
                        self.par['a_' + currSp +'_'  + '_'.join(list(combinationOfRegulators))] = \
                            int(boolodespace['boolval'])
                    elif self.settings['modeltype'] == 'heaviside':
                        self.par['w_' + currSp +'_'  + '_'.join(list(combinationOfRegulators))] = \
                            self.kineticParameterDefaults['heavisideOmega']*utils.heavisideThreshold(boolodespace['boolval'])                    

            # Close expressions
            if self.settings['modeltype'] == 'hill':
                num += ' )'
                den += ' )'
                f = '(' + num + '/' + den + ')'
            elif self.settings['modeltype'] == 'heaviside':
                # In the case of heaviside expressions, to prevent
                # numerical blowup, we trucate the magnitude of the
                # regulatory terms
                exponent += ')'
                maxexp = '10.' # '100'
                f = '(1./(1. + np.exp(np.sign('+exponent+')*min(' +maxexp +',abs(' + exponent+ ')))))'
            
            if currSp in self.proteinlist:
                Production =  f
                Degradation = 'p_' + currSp
                self.varspecs['p_' + currSp] = 'signalingtimescale*(y_max*' + Production \
                                           + '-' + Degradation + ')'
            else:
                
                Production = 'm_'+ currSp + '*' + f
                Degradation = 'l_x_'  + currSp + '*x_' + currSp
                self.varspecs['x_' + currSp] =  Production \
                                           + '-' + Degradation
                # Create the corresponding translated protein equation
                self.varspecs['p_' + currSp] = 'r_'+currSp+'*'+'x_' +currSp + '- l_p_'+currSp+'*'+'p_' + currSp
                
        ##########################################################                                       
            
        # Initialize variables between 0 and 1, Doesn't matter.
        xvals = [1. for _ in range(len(self.genelist))]
        pvals = [20. for _ in range(len(self.proteinlist))]    
        ics = {}
    
        for node, xv in zip(self.genelist, xvals):
            ics['x_' + node] = xv
            ics['p_' + node] = 0
        for node, pv in zip(self.proteinlist, pvals):
            ics['p_' + node] = pv
    
        self.ModelSpec['varspecs'] = self.varspecs
        self.ModelSpec['pars'] = self.par
        self.ModelSpec['ics'] = ics
        
        self.varmapper = {i:var for i,var in enumerate(self.ModelSpec['varspecs'].keys())}
        self.parmapper = {i:par for i,par in enumerate(self.ModelSpec['pars'].keys())}

    def writeModelToFile(self):
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
        self.path_to_ode_model = self.settings['outprefix'] / 'model.py'#os.path.dirname(os.path.realpath(__file__))    

        with open(self.path_to_ode_model,'w') as out:
            out.write('#####################################################\n')
            out.write('import numpy as np\n')
            out.write('# This file is created automatically\n')
            out.write('def Model(Y,t,pars):\n')
            out.write('    # Parameters\n')
            par_names = sorted(self.ModelSpec['pars'].keys())
            for i,p in enumerate(par_names):
                out.write('    ' + p + ' = pars[' + str(i) + ']\n')
            outstr = ''
            out.write('    # Variables\n')
            for i in range(len(self.varmapper.keys())):
                out.write('    ' + self.varmapper[i] + ' = Y[' + str(i) + ']\n')
                outstr += 'd' + self.varmapper[i] + ','
            for i in range(len(self.varmapper.keys())):
                vdef = self.ModelSpec['varspecs'][self.varmapper[i]]
                vdef = vdef.replace('^','**')
                out.write('    d' + self.varmapper[i] + ' = '+vdef+'\n')

            out.write('    dY = np.array([' + outstr+ '])\n')
            out.write('    return(dY)\n')
            out.write('#####################################################')

    
    def writeParametersToFile(self):
        """
        Writes dictionary of parameters to file

        :param ModelSpec: Model definition dictionary containing the keys 'ics', 'pars', and 'varspec' 
        :type ModelSpec: dict 
        :param outPrefix: Prefix to output folder name
        :type outPrefix: str
        :param outname: User defined name for parameters file. Default is parameters.txt
        :type outname: str (Optional)
        """
        with open(str(self.settings['outprefix']) + '/parameters.txt','w') as out:
            out.write('# Automatically generated by BoolODE\n')
            for k, v in self.ModelSpec['pars'].items():
                out.write(k+'\t'+str(v) + '\n')    
