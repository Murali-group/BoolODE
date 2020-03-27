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
from BoolODE import  utils
from BoolODE import simulatorCore as simulator 
from importlib.machinery import SourceFileLoader

class generateModel:
    def __init__(self, settings, parameterInputsDF, parameterSetDF, interactionStrengthDF) -> None:
        """
        initialize.
        """
        self.modelpath = settings['path']
        self.parameterInputsPath = settings['parameter_inputs_path']
        self.outprefix = settings['outprefix']
        self.add_dummy = settings['add_dummy']
        self.max_parents = settings['max_parents']
        self.sampleStd = settings['sample_std']
        self.withRules = list()
        self.withoutRules = list()
        self.allnodes = set()
        self.varspecs = dict()
        self.genelist = list()
        self.proteinlist = list()
        self.nodeTypeDF = pd.DataFrame()
        self.parameterInputsDF = parameterInputsDF
        self.parameterSetDF = parameterSetDF
        self.interactionStrengthDF = interactionStrengthDF
        self.df = pd.DataFrame()

        # Read the model definition
        # 1. populate self.df
        # 2. store genelist, withRules, withoutRules, allnodes
        # 3. initialize varspecs
        self.readBooleanRules()
        
        # Create model parameters and assign values
        self.getParameters()

        # Finally, create the model dictionary
        self.generateModelDict()
        
    def readBooleanRules(self):#, path, parameterInputsPath, outPrefix='',
        #                         add_dummy=False, max_parents=1):
        """
        Reads a rule file from path and stores the rules in a dictionary
        
        Parameters
        ----------
        :param path: Path to Boolean Rule file
        :type path: str
        :param parameterInputsPath: Path to file containing input parameters to the model
        :type parameterInputsPath: str
        :param add_dummy: Experimental feature, adds dummy targets to each gene in the model. Useful to change the density of the network
        :type add_dummy: bool
        :param max_parents: Experimental feature, specifies number of parent genes per dummy gene added
        :type max_parents: int
    
        :returns:
            - df: Dataframe containing rules
            - withoutrules: list of nodes in input file without rules
        """
        self.df = pd.read_csv(self.modelpath, sep='\t', engine='python')
        
        self.withRules = list(df['Gene'].values)

        ## Extract all nodes from the boolean rules
        ## Note that not all nodes might have a rule attached to them
        for ind,row in df.iterrows():
            rhs = row['Rule']
            rhs = rhs.replace('(',' ')
            rhs = rhs.replace(')',' ')
            tokens = rhs.split(' ')
            reg = [t for t in tokens if t not in ['not','and','or','']]
            self.allnodes.update(set(reg))
    
        self.withoutRules = list(allnodes.difference(set(withRules)))
        
        for n in withoutRules:
            if not Path(self.parameterInputsPath):
                print(n, "has no rule, adding self-activation.")
                self.df = self.df.append({'Gene':n,'Rule':n}, ignore_index=True)
            else:
                print("Treating %s as parameter" % {n})
                
        if add_dummy:
            self.addDummyGenes()
        # Remove
        #self.species = set(self.df['Gene'].values)
        
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
                    # Remove
                    # self.varspecs['x_' + node] = ''
                    # self.varspecs['p_' + node] = ''
                    self.genelist.append(node)
                else:
                    print(row['Type'],"not recognized. Treating as gene.")
                    self.varspecs.update({'x_' + node: '',
                                          'p_' + node: ''})
                    # Remove
                    # varspecs['x_' + node] = ''
                    # varspecs['p_' + node] = ''
                    self.genelist.append(node)
            # Remove
            # typeSpecified = set(self.nodeTypeDF['Node'].values)
            # for node in self.withRules:
            #     if node not in typeSpecified:
            #         self.genelist.append(sp)
        #print(self.df)
    
    def addDummyGenes(self):
        ## Add dummy genes
        ## determinisitcally add parents
        print('using max_parents=' + str(self.max_parents))
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
        self.df.to_csv(outPrefix + 'rules-with-added-genes.csv')        
        
    def getParameters(self):
        """
        Create dictionary of parameters and values. Assigns
        parameter values by evaluating the Boolean expression
        for each variable.
    
        :param DF: Table of values with two columns, 'Gene' specifies target, 'Rule' specifies Boolean function
        :type DF: pandas DataFrame
        :param identicalPars: Passed to utils.getSaneNval to set identical parameters
        :type identicalPars: bool
        :param samplePars: Sample kinetic parameters using a Gaussian distribution  centered around the default parameters
        :type samplePars: bool
        :param parameterInputsDF: Optional table that specifies input parameter values. Useful to specify experimental conditions.
        :type parameterInputsDF: pandas DataFrame
        :param parameterSetDF: Optional table that specifies predefined parameter set.
        :type parameterSetDF: pandas DataFrame
        :param interactionStrengthDF: Optional table that specifies interaction strengths. When not specified, default strength is set to 1.
        :type interactionStrengthDF: pandas DataFrame
        :returns:
            pars : dict
                Dictionary of parameters 
        """
    
        kineticParameterDefaults = utils.loadParameterValues()
        
        ## Derived parameters:
        ## The threshold is calculated as the half max value of the species
        x_max = kineticParameterDefaults['mRNATranscription']/kineticParameterDefaults['mRNADegradation']
        y_max = x_max*(kineticParameterDefaults['proteinTranslation']/kineticParameterDefaults['proteinDegradation'])
        
        hillThreshold = y_max/2
        heavisideOmega = 2./y_max
    
        parameterNamePrefixAndDefaultsAll = {
            # Hill coefficients
            'n_':kineticParameterDefaults['hillCoefficient'],
            # Thresholds
            'k_':kineticParameterDefaults['hillThreshold'],
            # Soft-Heaviside function steepness
            'sigmaH_':heavisideSigma
        }     
        parameterNamePrefixAndDefaultsGenes = {
            # mRNA transcription rates
            'm_':kineticParameterDefaults['mRNATranscription'],
            # mRNA degradation rates
            'l_x_':kineticParameterDefaults['mRNADegradation'],
            # Protein translation rate
            'r_':kineticParameterDefaults['proteinTranslation'],
            # protein degradation rates
            'l_p_':kineticParameterDefaults['proteinDegradation']
        }     
        
        par = dict()
        ## If there is even one signaling protein,
        ## create the y_max parameter
        if len(proteinlist) > 0:
            self.par['y_max'] = y_max
            self.par['signalingtimescale'] = signalingtimescale
            
        parameterInputsDict  = {}
        interactionStrengths = {}
    
        ## Carry out a series of checks to handle user defined input to model

        ## Check 1: handle "input parameters"
        if not self.parameterInputsDF.empty:
            # TODO 'parameter inputs' is confusing/misleading
            self.addParameterInputs(hillThreshold)

        ## Check 2: include user-specified interaction strengths
        if not self.interactionStrengthDF.empty:
            self.addInteractionStrengths(hillThreshold, interactionStrength)

        ## Assigining values to each kinetic parameter
        ## Create a parameter of each parameter type for each gene.
        ## If samplePars=True, sample values from normal distributions,
        ## else, set all parameters to defaults
        if self.samplePars:
            assignSampledParameterValues()
        else:
            assignDefaultParameterValues()

        # Final check:
        # Update user defined parameter values
        if not self.parameterSetDF.empty:
            for pname, pval in self.parameterSetDF.iterrows():
                self.par[pname] = pvalue[1]
                
    def addParameterInputs(self, hillThreshold):
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
        ## TODO Add appropriate parameters if heaviside            
        for i, row in self.parameterInputsDF.iterrows():
            # initialize
            parameterInputsDict[i] = {p:0 for p in withoutRules} ## Suspicious
            inputParams = ast.literal_eval(row['Inputs'])
            inputValues = ast.literal_eval(row['Values'])
            # Set to a max value
            parameterInputsDict[i].update({p:hillThreshold*2 for p,v in\
                                       zip(inputParams,inputValues) if v > 0})
        # Add input parameters, set to 0 by default
        self.par.update({p:0 for p in self.withoutRules})  ## Suspicious
        self.par.update({'k_'+p:hillThreshold for p in parameterInputsDict[0].keys()})
        self.par.update({'n_'+p:hillCoefficient for p in parameterInputsDict[0].keys()})


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
            self.par.update({'k_' + regulator + '_' + target: hillThreshold/interactionStrength})
        

    def assignDefaultParameterValues(self,
                                     parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes):
        print("Fixing rate parameters to defaults")
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            for node in self.withRules:
                self.par[parPrefix + node] = parDefault
                
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            for node in self.withRules:
                if node in genelist:
                    self.par[parPrefix + node] = parDefault


    def assignSampledParameterValues(self,
                                     parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes):
        """ 
        Sample kinetic parameters from a normal distribution, but in a range +- 10% of 
        the default value.
        """
        print("Sampling parameter values")
        print("Using std=" + str(self.sampleStd))
        lomult = 0.9
        himult = 1.1
        # Sample hillthreshold and hillcoeff
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsAll.items():
            sampledParameterValues = utils.getSaneNval(len(self.withRules),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=self.sampleStd*parDefault,\
                                                 identicalPars=identicalPars)
            for node, sparval in zip(self.withRules, sampledParameterValues):
                if node in self.genelist:
                    self.par[parPrefix + node] = sparval
                    
        transcriptionRate = 0.0 # Explain
        mRNADegradationRate = 0.0 # Explain
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            sampledParameterValues = utils.getSaneNval(len(species),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=sigmamult*parDefault,\
                                                 identicalPars=identicalPars)
            # This is very unclear! Why are we doing this?
            if identicalPars:
                if parPrefix == 'm_':
                    transcriptionRate = sampledParameterValues[0]
                if parPrefix == 'l_x_':
                    mRNADegradationRate = sampledParameterValues[0]
            else:
                if parPrefix == 'm_':
                    transcriptionRate = parDefault
                if parPrefix == 'l_x_':
                    mRNADegradationRate = parDefault
            
            for node, sparval in zip(self.withRules, sampledParameterValues):
                if node in genelist:
                    self.par[parPrefix + node] = sparval
        ## Based on these samples, estimate max value a
        ## gene can take at _steady-state_
        x_max = (transcriptionRate/mRNADegradationRate)         # What is this used for?

                    
    def generateModelDict(self):
        # , DF,identicalPars,
        #                   samplePars,
        #                   sampleStd,
        #                   withoutRules,
        #                   speciesTypeDF,
        #                   parameterInputsDF,
        #                   parameterSetDF,
        #                   interactionStrengthDF,
        #                   modeltype='hill'):
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
    
        if not self.parameterInputsDF.empty:
            inputs = set(self.withoutRules)
            for k in self.parameterInputs[0].keys():
                # Initialize variables to 0
                tempStr = k + " = 0"  
                exec(tempStr, boolodespace)
                
        for i,row in self.df.iterrows():
            # Basal expression:
            # Execute the rule to figure out
            # the value of alpha_0 or omega_0
            exec('booleval = ' + row['Rule'], boolodespace)         
            if modeltype == 'hill':
                par['alpha_'+row['Gene']] = int(boolodespace['booleval'])
            elif modeltype == 'heaviside':
                par['omega_' + row['Gene']] = utils.heavisideThreshold(boolodespace['booleval'])
        # End initialization

        # Assign values to logic parameters, and construct expressions
        for i,row in self.df.iterrows():
            ## Parse Boolean rule to get list of regulators
            allreg, regSpecies, regInputs = utils.getRegulatorsInRule(row['Rule'],
                                                                species,
                                                                inputs)

            # Basal expression term
            currSp = row['Gene']
            if modeltype == 'hill':
                num = '( alpha_' + currSp
                den = '( 1'
            elif modeltype == 'heaviside':
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
                    regulatorExpression = utils.createRegulatoryTerms(combinationOfRegulators,
                                                                      modeltype,
                                                                      strengthSpecified,
                                                                      regulatorsWithStrength,
                                                                      regSpecies)
                    if modeltype == 'hill':
                        # Create Numerator and Denominator
                        den += ' +' +  regulatorExpression
                        num += ' + a_' + currSp +'_'  + '_'.join(list(combinationOfRegulators)) + '*' + regulatorExpression
                    elif modeltype == 'heaviside':
                        exponent += ' + w_' + currSp + '_' + '_'.join(list(combinationOfRegulators)) +'*' + regulatorExpression
    
                    # evaluate rule to assign values to parameters
                    ##################################################
                    # for i1, row1 in DF.iterrows():                 #
                    #     exec(row1['Gene'] + ' = 0', boolodespace)  #
                    for node in withRules:                 #
                        exec(node + ' = 0', boolodespace)  #                        
                    # Set each regulator to ON, evaluate rule. we are looping
                    # over all such combinations of regulators
                    for geneInList in combinationOfRegulators:     #
                        exec(geneInList + ' = 1', boolodespace)    #
                                                                   #
                    exec('boolval = ' + row['Rule'], boolodespace) #
                    ##################################################
    
                    if modeltype == 'hill':
                        self.par['a_' + currSp +'_'  + '_'.join(list(combinationOfRegulators))] = \
                            int(boolodespace['boolval'])
                    elif modeltype == 'heaviside':
                        self.par['w_' + currSp +'_'  + '_'.join(list(combinationOfRegulators))] = \
                            heavisideOmega*utils.heavisideThreshold(boolodespace['boolval'])                    

            # Close expressions
            if modeltype == 'hill':
                num += ' )'
                den += ' )'
                f = '(' + num + '/' + den + ')'
            elif modeltype == 'heaviside':
                # In the case of heaviside expressions, to prevent
                # numerical blowup, we trucate the magnitude of the
                # regulatory terms
                exponent += ')'
                maxexp = '10.' # '100'
                f = '(1./(1. + np.exp(np.sign('+exponent+')*min(' +maxexp +',abs(' + exponent+ ')))))'
            
            if currSp in proteinlist:
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
        xvals = [1. for _ in range(len(genelist))]
        pvals = [20. for _ in range(len(proteinlist))]    
        ics = {}
    
        for node, xv in zip(self.genelist, xvals):
            ics['x_' + node] = xv
            ics['p_' + node] = 0
        for node, pv in zip(self.proteinlist, pvals):
            ics['p_' + node] = pv
    
        ModelSpec = {}
        ModelSpec['varspecs'] = varspecs
        ModelSpec['pars'] = par
        ModelSpec['ics'] = ics
        return ModelSpec, parameterInputs, genelist, proteinlist, x_max
    

    def getInitialCondition(self, ss, ModelSpec, rnaIndex,
                            proteinIndex,
                            genelist, proteinlist,
                            varmapper,revvarmapper):
        """
        Calculate the initial values of all state variables. 
        Takes into consideration user defined initial conditions, and computes the steady 
        states of the protein variables based on the estimated values of their corresponding genes.
    
        :param ss: Steady state array
        :type ss: ndarray
        :param ModelSpec: Dictionary of dictionary specifying the ODE model, containing parameters, initial conditions and equations.
        :type ModelSpec: dict
        :param rnaIndex: list of indices of genes
        :type rnaIndex: list
        :param proteinIndex: List of indices of proteins
        :type proteinIndex: list
        :param genelist: List of names of all genes in the model
        :type genelist: list
        :param proteinlist: List of names of all proteins in the model
        :type proteinlist: list
        :param varmapper: Mapper: {variable name : index}
        :type varmapper: dict
        :param revvarmapper: Mapper: {index : variable name}
        :type revvarmapper: dict
        :returns:
            - newics: List containing new initial conditions
        """
        
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
    
