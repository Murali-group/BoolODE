import pandas as pd
from tqdm import tqdm
import numpy as np
import os
import sys
from pathlib import Path
from optparse import OptionParser

def parseArgs(args):
    parser = OptionParser()
    
    parser.add_option('', '--outPrefix', type = 'str',default='',
                      help='Prefix for output files. This defaults to the value of --input-path.')
    
    parser.add_option('-p', '--input-path', type='str',
                      help='Path to simulation output folder')
    
    parser.add_option('-n', '--nCells', type='int',default='100',
                      help='Number of cells in each dataset. This should be less than the  number of files under /simulations')
    
    parser.add_option('-d', '--nDatasets', type='int',default='1',
                      help='Number of datasets to be generated')        
    
    (opts, args) = parser.parse_args(args)

    return opts, args

def genSamples(opts):
    if opts.input_path is None:
        print('Input path is invalid.')
        sys.exit()
    clusterdf = pd.read_csv(opts.input_path + '/ClusterIds.csv', index_col=0)
    numclusters = clusterdf.cl.unique()
    num_experiments = clusterdf.shape[0]
    if opts.nCells > num_experiments:
        print('nCells should be less than num of experiments')
        sys.exit()
    df = pd.read_csv(opts.input_path + '/simulations/E0.csv',index_col=0)
    maxtime = len(df.columns)
    for did in range(1, opts.nDatasets + 1):
        # example:
        # Beeline/inputs/DYN-LI-500-1/...
        outfpath = opts.input_path + '/' + opts.outPrefix + '-' + str(opts.nCells) + '-' + str(did) 
        if not os.path.exists(outfpath):
            print(outfpath, "does not exist, creating it...")
            os.makedirs(outfpath)
        # Create cell ids
        simids = np.random.choice(range(num_experiments), size=opts.nCells, replace=False)
        fids = ['E'+ str(sid) + '.csv' for sid in simids]            
        timepoints = np.random.choice(range(1,max_time), size=opts.nCells)
        min_t = min(timepoints)
        max_t = max(timepoints)
        pts = [(t - min_t)/(max_t - min_t) for t in timepoints]
        cellids = ['E' + str(sid) + '_' + str(t) for sid, t in zip(simids, timepoints)] 
        # Read simulations from input dataset #psetid
        # to build a sample
        sample = []
        for fid, cid in tqdm(zip(fids, cellids)):
            df = pd.read_csv( opts.input_path + '/simulations/' + fid, index_col=0)
            df.sort_index(inplace=True)
            sample.append(df[cid].to_frame())
        sampledf = pd.concat(sample,axis=1)
        sampledf.to_csv(outfpath + 'ExpressionData.csv')
        ## Read refNetwork.csv
        refdf = pd.read_csv(opts.input_path + '/refNetwork.csv')
        refdf.to_csv(outfpath + 'refNetwork.csv',index=False)
        if len(numclusters) == 1:
            ptdf = pd.DataFrame(np.array(pts),
                            index=pd.Index(cellids),columns = ['PseudoTime'])
        else:
            subsetcluster = clusterdf.loc[[cid.split('_')[0] for cid in cellids]]
            ptdf = pd.DataFrame(index=pd.Index(cellids),
                                columns = ['PseudoTime' + str(1+i) for i in range(numclusters)])
            for (cid, row), pt in zip(ptdf.iterrows(), pts):
                ptdf.loc[cid]['PseudoTime' +\
                              str(int(subsetcluster.loc[cid.split('_')[0]]['cl']) + 1)] = round(pt,5)
        ptdf.to_csv(outfpath + 'PseudoTime.csv',na_rep='NA')                    

def main(args):
    opts, args = parseArgs(args)
    genSamples(opts)

if __name__ == '__main__':
    main(sys.argv)
    
