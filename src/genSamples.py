import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from optparse import OptionParser 
import sys
seed = 0
np.random.seed(seed)

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
    
    
    parser.add_option('-s', '--nSamples', type='int',default='10',
                      help='Number of random samples of size n.')    
    
    (opts, args) = parser.parse_args(args)

    return opts, args

def genSamples(opts):
    ExprDF = pd.read_csv(opts.expr,index_col=0, header = 0)
    ptDF = pd.read_csv(opts.pseudo,index_col=0, header = 0)    
    refDF = pd.read_csv(opts.refNet, index_col = None, header = 0)
    
    # figure out which experiments to drop
    headers = ExprDF.columns
    experiments = set([h.split('_')[0] for h in headers])
    dfmax = ExprDF.values.max()
    droplist = []
    for e in experiments:
        ecount = 0
        cellcount = 0
        for col in  ExprDF.columns:
            if e == col.split('_')[0]:
                cellcount += 1
                if ExprDF[col].max() < 0.1*dfmax:
                    ecount +=1
        if ecount > 0:
            droplist.append(e)

    for d in tqdm(droplist):
        for col in ExprDF.columns:
            if d == col.split('_')[0]:
                ExprDF.drop(col,axis=1,inplace=True)

    for i in range(opts.nSamples):
        path = opts.outPrefix + '_' +str(opts.nCells) +'_' + str(i)
        if not os.path.exists(path):
            os.makedirs(path)

        refDF.to_csv(path + '/refNetwork.csv',index=False)

        # New samples
        SampleDF = ExprDF.sample(n = opts.nCells, axis = 'columns')
        SampleDF.to_csv(path +'/ExpressionData.csv')
        samplePtDF = ptDF.loc[SampleDF.columns,:]
        samplePtDF.to_csv(path +'/PseudoTime.csv')
        
def main(args):
    opts, args = parseArgs(args)
    genSamples(opts)
    
if __name__ == "__main__":
    main(sys.argv)
