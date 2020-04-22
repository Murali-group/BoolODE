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
    """Class that holds model attributes. Provides helper functions to convert a Boolean model
    to an ODE model
    
    :param settings: User defined job settings. This is created in boolode.do_post_processing()
    :type settings: dict
    :param parameterInputsDF: User specified parmeter inputs. Specified in config YAML as a path.
    :type parameterInputsDF: pandas DataFrame
    :param parameterSetDF: User specified kinetic parameters. Specified in config YAML as a path.
    :type parameterSetDF: pandas DataFrame
    :param interactionStrengthDF: User specified interaction strengths. Specified in config YAML as a path
    :type interactionStrengthDF: pandas DataFrame
    """
    def __init__(self, settings, parameterInputsDF, parameterSetDF, interactionStrengthDF) -> None:
        """
        This class is initialized with the job_settings dictionary created 
        by the BoolODE class. Additionally, three dataframe objects are required:

        1. parameterInputsDF - specifies parameter inputs
        2. parameterSetDF - specifies pregenerated python file containing kinetic parameter values
        3. interactionStrengthDF - specifies interactions strengths

        By default, if these aren't specified by the user, BoolODE stores empty pandas DataFrame objects.
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
        Reads a rule file from path and stores the rules in a dictionary.
        Performs the following transformations:

        1. Parses each Boolean rule to extract nodes in the model
        2. Identifies nodes without corresponding rules: (a) If a parameterInput file is specified and the node under consideration
           is present in the file, the node is treated as an input parameter
           (b) Else, a self-edge is added to the node
        3. Adds dummy genes - *This is experimental*
        4. Finally, create a gene (x\_) and protein (p\_)  variable corresponding to 
           each node in the network, taking into consideration user defined types.
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
        """
        Add dummy genes. 
        Starts with the user specified Boolean model, and grows the
        network by adding new 'dummy' nodes to each existing node.
        In the current configuration, the user can specify `max_parents`
        which specifies the number of parent nodes of every dummy node.
        Currently, twice as many dummy nodes as nodes in the original graph
        are added. 

        .. todo::
            Expand this function by adding user defined number of nodes.

        .. warning::
            This feature is experimental.
        """
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
        """Create and assigns the kinetic parameters for each equations the
        corresponding parameter value.

        .. note::
            Default kinetic parameter values are stored in ./parameters.yaml. This
            is provided to make it easy for users to define different defaults, but it 
            other parameter values have not been tested extensively, and might not
            yield the expected results.
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
            # Soft-Heaviside function steepness,
            # This is to allow for flexibility in the future.
            # not currently used
            # 'sigmaH_':self.kineticParameterDefaults['heavisideSigma']
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

        1. 'Input' specifies the name of input parameter
        2. 'Value' specified the input state, typically 1 or 0, but for HIGH and LOW inputs respectively, but any float between 0 and 1 is valid input.
          
        Now, we create a fake hill term for this input, even
        though there is no equation corresponding to it. Thus,
        we create a hillThreshold and hillCoefficient 
        term for it. This is useful so we can treat this input
        parameter as a "regulator" and leave the structure of
        the logic term unchanged.

        .. todo::
            Add appropriate parameters if heaviside
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
        """Check if interaction strengths are specified.
        Create a new parameter specifying the strength
        of a given interaction. For all other equations in
        which "regulator" appears in, its hillThreshold value
        will remain the default value.

        .. todo::
            Add appropriate parameters if heaviside
        """
        for i,row in self.interactionStrengthDF.iterrows():
            regulator = row['Gene2']
            target = row['Gene1']
            interactionStrength = row['Strength']
            if self.settings['modeltype'] == 'hill':
                self.par.update({'k_' + regulator + '_' + target: hillThreshold/interactionStrength})

    def assignDefaultParameterValues(self,
                                     parameterNamePrefixAndDefaultsAll,
                                     parameterNamePrefixAndDefaultsGenes):
        """
        Set each kinetic parameter to its default value.
        """
        print("Fixing rate parameters to defaults")
        if self.settings['modeltype'] == 'hill':
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
        Sample kinetic parameters from a truncated normal distribution with mean=default
        parameter value, and standard deviation = sample_std, in a range +- 10% of 
        the default value. 
        """
        print("Sampling parameter values")
        print("Using std=" + str(self.settings['sample_std']))
        lomult = 0.9
        himult = 1.1
        if self.settings['modeltype'] == 'hill':
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
                    
        for parPrefix, parDefault in parameterNamePrefixAndDefaultsGenes.items():
            sampledParameterValues = utils.getSaneNval(len(self.withRules),\
                                                 lo=lomult*parDefault,\
                                                 hi=himult*parDefault,\
                                                 mu=parDefault,\
                                                 sig=self.settings['sample_std']*parDefault,\
                                                 identicalPars=self.settings['identical_pars'])
            
            for node, sparval in zip(self.withRules, sampledParameterValues):
                if node in self.genelist:
                    self.par[parPrefix + node] = sparval

    def createRegulatoryTerms(self, currgene, combinationOfRegulators,
                              regSpecies):
        """Creates the terms in the activation function which governs
        the transcriptional state of the gene. 
        If modeltype is 'hill', each term is composed of products of Hill functions, ranging over
        all possible combinations of regulators. 
        Similarly, if modeltype is 'heaviside', each term is simply a
        product of the activities of the combinations of regulators, for all combinations.
        Here, the form of a single term is generated for a single combination of regulators.
        
        :param currgene: Name of the current gene
        :type currgene: str
        :param combinationOfRegulators: a list of all combinations of regulators of currgene
        :type combinationOfRegulators: list
        """
        strengthSpecified = False
        if not self.interactionStrengthDF.empty:
            if currgene in self.interactionStrengthDF['Gene1'].values:
                strengthSpecified = True
                regulatorsWithStrength = set(self.interactionStrengthDF[\
                                                                        self.interactionStrengthDF['Gene1']\
                                                                        == currgene]['Gene2'].values)            
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

    def shsparse(self, regulators, rule):
        parser = BoolParser(str(rule))
        shsexp = parser.constructPolynomial(parser.createBoolTree())
        #basal = '-0.3 + ' # This works, but I don't know which threshold will approximate the Hill function behavior better
        # sigma = 10        
        basal = '-0.5 + '
        sigma = 10
        ymax = self.kineticParameterDefaults['y_max']
        # NOTE: Another possibility here is to raise each term to an exponent
        # When I tried this out, most simulations "failed" by producing close-to-zero expression for all genes
        # This might yet work in some regime of (sigma, exponent, noise strength) which I have not yet come across
        for v in regulators:
            shsexp = shsexp.replace(v, \
                                    '(p_' + v + '/' + str(ymax) + ')') #'((p_' + v + '/' + str(ymax) + ')**3.0)')
        shsexp = basal + shsexp
        shsexp = '(1.0/(1 + np.exp(-'+str(sigma)+'*(' + shsexp +'))))'
        return shsexp
    
        
    def generateModelDict(self):
        """
        Take a DataFrame object with Boolean rules,
        construct ODE equations for each variable.
        This is the core function in BoolODE.
        In a nutshell, a temporary namespace called boolodespace
        is created. In this namespace, the nodes  in the Boolean network
        are first initialized to 0, or OFF. Then, for each gene/node,
        its corresponding rule is evaluated by setting its regulators to
        all combinations binary states. The outcome of each of these Boolean
        rule evaluations is used to decide the value of the activation strength
        parameter 'a' in Hill functions. For the Heaviside functions, the boolean rule
        does not need to be evaluated. Rather, the rule is expressed directly
        into a polynomial, by the BoolParser() class.
        Each step of this conversion is documented in the code directly.
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
        # End initialization

        # Assign values to logic parameters, and construct expressions
        for i,row in self.df.iterrows():
            ## Parse Boolean rule to get list of regulators
            allreg, regSpecies, regInputs = utils.getRegulatorsInRule(row['Rule'],
                                                                      self.withRules,
                                                                      self.inputs)

            # Basal expression term
            currgene = row['Gene']
            if self.settings['modeltype'] == 'hill':
                num = '( alpha_' + currgene
                den = '( 1'
                # Loop over combinations of regulators        
                for i in range(1,len(allreg) + 1):
                    for combinationOfRegulators in combinations(allreg,i):
                        regulatorExpression = self.createRegulatoryTerms(currgene, combinationOfRegulators,
                                                                          regSpecies)
                        # Create Numerator and Denominator
                        den += ' +' +  regulatorExpression
                        num += ' + a_' + currgene +'_'  + '_'.join(list(combinationOfRegulators)) + '*' + regulatorExpression
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
                        self.par['a_' + currgene +'_'  + '_'.join(list(combinationOfRegulators))] = \
                            int(boolodespace['boolval'])
                # Close expressions
                num += ' )'
                den += ' )'
                f = '(' + num + '/' + den + ')'
                
            elif self.settings['modeltype'] == 'heaviside':
                f = self.shsparse(allreg, row['Rule'])
            
            if currgene in self.proteinlist:
                Production =  f
                Degradation = 'p_' + currgene
                self.varspecs['p_' + currgene] = 'signalingtimescale*(y_max*' + Production \
                                           + '-' + Degradation + ')'
            else:
                Production = 'm_'+ currgene + '*' + f
                Degradation = 'l_x_'  + currgene + '*x_' + currgene
                self.varspecs['x_' + currgene] =  Production \
                                           + '-' + Degradation
                # Create the corresponding translated protein equation
                self.varspecs['p_' + currgene] = 'r_'+currgene+'*'+'x_' +currgene + '- l_p_'+currgene+'*'+'p_' + currgene
                
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
        Writes model to file as a python function.
        The ODE model generated using generateModelDict() is defined as 
        a python ODE function called Model(). Model() takes 3 arguments:

        1. The current model state vector Y
        2. The current time t
        3. A list of parameters pars

        The function returns a vector of time derivatives computed from the ODEs.
        Model() is written model.py in the directory of the current job
        """
        self.path_to_ode_model = self.settings['outprefix'] / 'model.py'

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
        Writes dictionary of parameters to file. 
        """
        with open(str(self.settings['outprefix']) + '/parameters.txt','w') as out:
            out.write('# Automatically generated by BoolODE\n')
            for k, v in self.ModelSpec['pars'].items():
                out.write(k+'\t'+str(v) + '\n')    


class Node():
    def __init__(self, left, op, right, leftstr='', rightstr=''):
        self.left = left
        self.right = right
        self.op = op
        self.rightstr = rightstr
        self.leftstr = leftstr

class BoolParser():
    def __init__(self, expr):
        self.expr = expr
        self.operators = ['and', 'or', 'not']
        self.tokenlist = self.preprocessExpression()

    def preprocessExpression(self):
        s = self.expr
        s = s.replace('(', ' ( ')
        s = s.replace(')', ' ) ')
        # remove trailing whitespace
        s = s.strip()
        tokenlist = [t for t in s.split(' ') if t != '']
        return tokenlist

    def createBoolTree(self, tokenlist=None):
        if tokenlist is None:
            tokenlist = self.tokenlist
        tokenlist = self.removeDelimiters(tokenlist)
        delimiters = self.getDelimiterPositions(tokenlist)
        tokenind = 0
        while tokenind < len(tokenlist):
            # print(tokenlist, tokenind)
            if len(delimiters) > 0:
                for o, c in delimiters:
                    if o == 0:
                        break
                if o == 0:
                    tokenind = c + 1
            if tokenlist[tokenind] in self.operators:
                break
            else:
                tokenind += 1
        if tokenind == len(tokenlist):
            return Node('', tokenlist[-1], '')
        t = tokenlist[tokenind]
        if t == 'not':
            left = ''
            right = tokenlist[tokenind+1:]
        else:
            left = tokenlist[:tokenind]
            right = tokenlist[tokenind+1:]
        
        argumentdict = {}
        for k, v in zip(['left','right'],[left,right]):
            if type(v) is str:
                argumentdict[k] = v
                argumentdict[k+'str'] = ''
            elif type(v) is list and len(v) == 1:
                argumentdict[k] = v[0]
                argumentdict[k+'str'] = ''        
            else:
                argumentdict[k+'str'] =  '(' + ' '.join(v) + ')'
                argumentdict[k] = self.createBoolTree(v)
        return Node(argumentdict['left'], t, argumentdict['right'],
                    leftstr=argumentdict['leftstr'],
                    rightstr=argumentdict['rightstr'])

    def getParenCount(self, tokenlist):
        parencount = 0
        for t in tokenlist:
            if t == '(':
                parencount += 1
            if t == ')':
                parencount -= 1
        return parencount

    def removeDelimiters(self, tokenlist):
        delimiters = self.getDelimiterPositions(tokenlist)
        if len(delimiters) > 0:
            for o, c in delimiters:
                if o == 0:
                    break
            if o == 0 and len(tokenlist) == c + 1 + o:
                tokenlist = list(tokenlist[1:-1])
                delimiters = self.getDelimiterPositions(tokenlist)
                # Call itself again
                tokenlist = self.removeDelimiters(list(tokenlist))
        return(tokenlist)

    def printBoolTree(self, currnode):
        l = currnode.left
        r = currnode.right

        if type(currnode.left) is not str:
            l = currnode.leftstr
            self.printBoolTree(currnode.left)

        if type(currnode.right) is not str:
            r = currnode.rightstr
            self.printBoolTree(currnode.right)
        print(l, currnode.op, r)

    def constructPolynomial(self, currnode):
        l = currnode.left
        r = currnode.right

        if type(currnode.left) is not str:
            l = self.constructPolynomial(currnode.left)

        if type(currnode.right) is not str:
            r= self.constructPolynomial(currnode.right)
        expr = ''
        if currnode.op == 'or':
            expr = '(1. - (1. - ' +l + ')*(1. - ' + r + '))'
        elif currnode.op == 'and':
            expr = '(' + l +'*' + r + ')'
        elif currnode.op == 'not':
            expr = '(1. - ' + r + ')'
        else:
            expr = currnode.op
        return expr

    def getDelimiterPositions(self, tokenlist):
        separatordict = {'open':[],'close':[]}
        for i, t in enumerate(tokenlist):
            if t == '(':
                separatordict['open'].append(i)
            elif t == ')':
                separatordict['close'].append(i)
        if len(separatordict['open']) != len(separatordict['close']):
            print(separatordict, tokenlist)
            print('Imbalanced expession!')
            sys.exit()
        explocations = []
        separatordict['close'].reverse()
        while len(separatordict['open']) > 0:
            c = separatordict['close'].pop()
            o = max([l for l in separatordict['open'] if l < c])
            separatordict['open'].remove(o)
            explocations.append((o,c))
        
        return explocations

