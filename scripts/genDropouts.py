import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from optparse import OptionParser
import time

def parseArgs(args):
    parser = OptionParser()
    
    parser.add_option('', '--outPrefix', type = 'str',default='',
                      help='Prefix for output files.')
    
    parser.add_option('-e', '--expr', type='str',
                      help='Path to expression data file')
    
    parser.add_option('-p', '--pseudo', type='str',
                      help='Path to pseudotime file')
        
    parser.add_option('-r', '--refNet', type='str',
                      help='Path to reference network file')
    
    parser.add_option('-n', '--nCells', type='int',default='100',
                      help='Number of cells to sample.')    
    
    parser.add_option('-d', '--dropout', action='store_true',default=False,
                      help='Carry out dropout analysis?')
    
    parser.add_option('', '--drop-cutoff', type="float",default=0.5,
                      help='Specify quantile cutoff on gene expression')    
    
    parser.add_option('', '--drop-prob', type='float',default=0.5,
                      help='Specify the probability of dropping a gene below quantile q ')    
    
    parser.add_option('-i', '--samplenum', type='int',default=None,
                      help='Sample Number')
            
    (opts, args) = parser.parse_args(args)

    return opts, args

def genSamples(opts):
    if opts.samplenum is None:
        print('Please specify sample number')
        sys.exit()
    if opts.dropout:
        dropoutCutoffs = opts.drop_cutoff
    else:
        dropoutCutoffs = 0

    ## Read the ExpressionData.csv file
    expDF = pd.read_csv(opts.expr, index_col=0)
    ## Read the PseudoTime.csv file
    PTDF = pd.read_csv(opts.pseudo, index_col=0)
    ## Read the refNetwork.csv file
    refDF = pd.read_csv(opts.refNet, index_col=0)

    ## Generate output path for the dropout datasets
    path = opts.outPrefix + '-' +str(opts.nCells) +'-' + str(opts.samplenum)
    if dropoutCutoffs != 0:
        path += '-' + str(int(100*dropoutCutoffs))+  '-' + str( str(opts.drop_prob))
    if not os.path.exists(path):
        os.makedirs(path)
    # Dropout here
    # copy over PT and refNetwork files
    refDF.to_csv(path + '/refNetwork.csv')
    PTDF.to_csv(path+'/PseudoTime.csv')        
    # Drop-out genes if they are less than the 
    # percentile value @ "dc" with 50% chance
    DropOutDF = expDF.copy()
    if dropoutCutoffs != 0:
        quantileExp = expDF.quantile(q = dropoutCutoffs, axis = 'columns')
        for idx, row in tqdm(expDF.iterrows()):
            for col in expDF.columns:
                if row[col] < quantileExp.loc[idx]:
                    cointoss = np.random.random()
                    if cointoss < opts.drop_prob:
                        DropOutDF.loc[idx,col] = 0.0

    DropOutDF.to_csv(path + '/ExpressionData.csv')
        
def main(args):
    opts, args = parseArgs(args)
    genSamples(opts)

if __name__ == '__main__':
    main(sys.argv)
