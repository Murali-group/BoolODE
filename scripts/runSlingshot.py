__author__ = 'Aditya Pratapa'
"""
Note: This script runs Slingshot in order to compute
the pseudotime for a set of cells simulated using BoolODE.
In order for this script to work, please ensure that the 
Slingshot docker is working correctly. 
"""
import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from optparse import OptionParser 
import seaborn as sns

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
    
    parser.add_option('-n', '--nCells', type='int',default='100',
                      help='Number of cells to sample.')    
    
    parser.add_option('-s', '--nSamples', type='int',default='10',
                      help='Number of random samples of size n.')
    
    parser.add_option('-d', '--dropout', action='store_true',default=False,
                      help='Carry out dropout analysis?')
            
    parser.add_option('-c', '--nClusters', type='int',default='1',
                      help='Number of expected clusters in the dataset.')    
    
    parser.add_option('', '--noEnd', action='store_true',default= False,
                      help='Do not force SlingShot to have an end state.')    
    
        
    parser.add_option('-r', '--perplexity', type='int',default=500,
                      help='Perplexity for tSNE.')
    
    
    (opts, args) = parser.parse_args(args)

    return opts, args

        

def computeSSPT(ExpDF, ptDF, nClust, outPaths, noEnd = False, perplexity = 500):
    '''
    Compute PseudoTime using 'slingshot'.
    Needs the input GenexCells expression data frame.
    Needs number of clusters to be expected in the 
    dataset. 
    E.g., Linear: k=1 (no PseudoTime inference is done)
    E.g., Bifurcating: k=3 (1 initial and 2 terminal)
    E.g., Trifurcating: k=4 (1 initial and 3 terminal)
    '''
    if nClust == 1:
        # Return simulation time as PseduoTime
        for outPath in outPaths:
            ptDF.loc[ExpDF.columns].to_csv(outPath+"/PseudoTime.csv", columns =['Time'])
 
    else:
        ### Compute PseudoTime ordering using slingshot
        
        # Step-1: Compute dimensionality reduction
        # Currently only does TSNE
        # TODO: Add PCA
        print("Computing TSNE...")
        DimRedRes = TSNE(n_components = 2, perplexity = perplexity).fit_transform(ExpDF.T)

        # Step-2: Convert TSNE results to a dataframe
        DimRedDF = pd.DataFrame(DimRedRes,columns=['dim1','dim2'],
                                 index=pd.Index(list(ExpDF.columns)))
        
        # Step-3: Compute kMeans clustering
        DimRedDF.loc[:,'cl'] = KMeans(n_clusters = nClust).fit(ExpDF.T).labels_
        # Identify starting cluster
        DimRedDF.loc[:,'pt'] = ptDF.min(axis='columns')
        
        # Step-4: Identify cells corresponding to initial and final states
        # Cells in initial states are identified from cluster with smallest
        # mean experimental time
        # Cells in final states are rest of the clusters
        startClust = DimRedDF.groupby('cl').mean()['pt'].idxmin()
        endClust = ','.join([str(ix) for ix in DimRedDF.groupby('cl').mean()['pt'].index if ix != startClust])
        startClust = str(startClust)
        
        # Step-5: Create a temporary directory and write files necessary to run slingshot
        os.makedirs(outPaths, exist_ok = True)

        DimRedDF.to_csv(outPaths + '/rd.tsv', columns = ['dim1','dim2'],sep='\t')
        DimRedDF.to_csv(outPaths + '/cl.tsv', columns = ['cl'],sep='\t')
        if noEnd:
            cmdToRun= " ".join(["docker run --rm -v", str(Path.cwd()) + "/" +outPaths +"/:/data/temp",
                    "slingshot:base /bin/sh -c \"Rscript data/run_slingshot.R",
                    "--input=/data/temp/rd.tsv --input-type=matrix",
                    "--cluster-labels=/data/temp/cl.tsv",
                    "--start-clus="+startClust+'\"'])

        else:
            cmdToRun= " ".join(["docker run --rm -v", str(Path.cwd())+ "/" +outPaths +"/:/data/temp",
                                "slingshot:base /bin/sh -c \"Rscript data/run_slingshot.R",
                                "--input=/data/temp/rd.tsv --input-type=matrix",
                                "--cluster-labels=/data/temp/cl.tsv",
                                "--start-clus="+startClust, "--end-clus="+endClust+'\"'])
        print(cmdToRun)
        os.system(cmdToRun)

        # os.system("cp temp/PseudoTime.csv "+outPaths+"/SlingshotPT.csv")
        # os.system("cp temp/curves.csv "+outPaths+"/curves.csv")
        # os.system("cp temp/rd.tsv "+outPaths+"/rd.tsv")
        # os.system("cp temp/cl.tsv "+outPaths+"/cl.tsv")
        
        # Do this only for the first file
        tn = pd.read_csv(outPaths+"/rd.tsv",
             header = 0, index_col =None, sep='\t')

        cl = pd.read_csv(outPaths+"/cl.tsv",
                         header = 0, index_col =0, sep='\t')

        curveFile = open(outPaths+"/curves.csv","r")
        curveLst = []
        lneCnt = 0
        for line in curveFile:
            curveLst.append([float(p) for p in line.strip().split(',')])
            lneCnt += 1

        tn['cl'] = cl.values
        tn.columns = ['CellID','dim1','dim2','kMeans']
        tn.index = tn.CellID

        f, axes = plt.subplots(2, 2, figsize=(7.5, 7.5))

        # Plot slingshot pseudotime 
        # and original clusters
        detPT = pd.read_csv(outPaths+"/SlingshotPT.csv",
             header = 0, index_col = 0)
        print()
        colNames = detPT.columns
        for colName in colNames:
            # Select cells belonging to each pseudotime trajectory
            index = detPT[colName].index[detPT[colName].notnull()]
            tn.loc[index,colName] = detPT.loc[index,colName]


            sns.scatterplot(x='dim1',y='dim2', 
                            data = tn.loc[index],  hue = colName,
                            palette = "viridis", 
                            ax = axes[1][0])
            plt.legend([])

        for line in range(0, lneCnt, 2):
            sns.lineplot(x= curveLst[line+1],y=curveLst[line],
                            color = "k", ax = axes[1][0])

        sns.scatterplot(x='dim1',y='dim2', 
                        data = tn,  hue = 'kMeans',
                        palette = "Set1", 
                        ax = axes[1][1])

        # Plot deterministic pseduotime 
        # and original clusters
        detPT = pd.read_csv(outPaths+"/PseudoTime.csv",
             header = 0, index_col = 0)
        colNames = detPT.columns
        for idx in range(len(colNames)):
            # Select cells belonging to each pseudotime trajectory
            colName = colNames[idx]
            index = detPT[colName].index[detPT[colName].notnull()]
            tn.loc[index,'Original'] = int(idx)

        tn['ExpTime'] = detPT.min(axis='columns')

        sns.scatterplot(x='dim1',y='dim2', 
                        data = tn,  hue = 'ExpTime',
                        palette = "viridis", 
                        ax = axes[0][0])

        #tn['Original'] = tn['Original'].astype('category')
        sns.scatterplot(x='dim1',y='dim2', 
        data = tn,  hue = 'Original',
        palette = "Set1", 
        ax = axes[0][1])

        axes[0][0].get_legend().remove()
        axes[0][0].title.set_text('Experiment Time')
        axes[0][1].get_legend().remove()
        axes[0][1].title.set_text('Original Trajectories')

        axes[1][0].get_legend().remove()
        axes[1][0].title.set_text('Slingshot Pseudotime')

        axes[1][1].get_legend().remove()
        axes[1][1].title.set_text('kMeans Clustering')

        f.tight_layout()

        tn.to_csv(outPaths+"/Updated_rd.tsv",
                  sep='\t')
        plt.savefig(outPaths+"/SlingshotOutput.png")
    os.system("rm -rf temp/")

def main(args):
    opts, args = parseArgs(args)
    
    ExprDF = pd.read_csv(opts.expr,index_col=0, header = 0)
    ptDF = pd.read_csv(opts.pseudo,index_col=0, header = 0)    

    # Compute PseudoTime using slingshot
    # TODO: Add other methods
    computeSSPT(ExprDF, ptDF, opts.nClusters, opts.outPrefix, opts.noEnd, opts.perplexity)
        
if __name__ == "__main__":
    main(sys.argv)
