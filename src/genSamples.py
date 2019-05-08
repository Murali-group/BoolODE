
import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE

from optparse import OptionParser 


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
    
        
    parser.add_option('-c', '--nClusters', type='int',default='1',
                      help='Number of expected clusters in the dataset.')    
    
            
    parser.add_option('', '--noEnd', action='store_true',default= False,
                      help='Do not force SlingShot to have an end state.')    
    
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
                
    dropoutCutoffs = [0,0.1,0.2,0.5]
    
    for i in range(opts.nSamples):
        # New samples
        SampleDF = ExprDF.sample(n = opts.nCells, axis = 'columns')
        # Dropout here

        meanExpression = {}
        for gene in SampleDF.index:
            meanExpression[gene] = SampleDF.loc[gene].mean()
            
        maxExp = max([v for g,v in meanExpression.items()])
        outPaths = []
        for dc in dropoutCutoffs:
            path = opts.outPrefix + '_' +str(opts.nCells) +'_' + str(i) + '_' + str(dc)
            if not os.path.exists(path):
                os.makedirs(path)
            outPaths.append(path)

        # Compute PseudoTime using slingshot
        # TODO: Add other methods
        computeSSPT(SampleDF, ptDF, opts.nClusters, outPaths, opts.noEnd)
        
        for dc,path in zip(dropoutCutoff,outPaths):
            path = opts.outPrefix + '_' +str(opts.nCells) +'_' + str(i) + '_' + str(dc)
            if not os.path.exists(path):
                os.makedirs(path)
                
            refDF.to_csv(path + '/refNetwork.csv',index=False)

            todrop = []
            for gene, meanGeneExp in meanExpression.items():
                if meanGeneExp <= dc*maxExp:
                    todrop.append(gene)
            DropOutDF = SampleDF.copy()
            for gene, row in SampleDF.iterrows():
                if gene in todrop:
                    for cell in DropOutDF.columns:
                        cointoss = np.random.random()
                        if cointoss < 0.5:
                            DropOutDF[cell].loc[gene] = 0.0
            DropOutDF.to_csv(path + '/ExpressionData-dropout-'+str(dc)+'.csv')

        #SampleDF.to_csv(path +'/ExpressionData.csv')
        

def computeSSPT(ExpDF, ptDF, nClust, outPaths, noEnd = False):
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
            ptDF.loc[ExpDF.columns,:].to_csv(outPath+"/PseduoTime.csv")
    else:
        ### Compute PseudoTime ordering using slingshot
        
        # Step-1: Compute dimensionality reduction
        # Currently only does TSNE
        # TODO: Add PCA
        DimRedRes = TSNE(n_components = 2).fit_transform(ExpDF.T)

        # Step-2: Convert TSNE results to a dataframe
        DimRedDF = pd.DataFrame(DimRedRes,columns=['dim1','dim2'],
                                 index=pd.Index(list(ExpDF.columns)))
        
        # Step-3: Compute kMeans clustering
        DimRedDF.loc[:,'cl'] = KMeans(n_clusters = nClust).fit(ExpDF.T).labels_
        DimRedDF.loc[:,'pt'] = ptDF.loc[DimRedDF.index,'Time']
        
        # Step-4: Identify cells corresponding to initial and final states
        # Cells in initial states are identified from cluster with smallest
        # mean experimental time
        # Cells in final states are rest of the clusters
        startClust = DimRedDF.groupby('cl').mean()['pt'].idxmin()
        endClust = ','.join([str(ix) for ix in DimRedDF.groupby('cl').mean()['pt'].index if ix != startClust])
        startClust = str(startClust)
        
        # Step-5: Create a temporary directory and write files necessary to run slingshot
        os.makedirs("temp/", exist_ok = True)

        DimRedDF.to_csv('temp/rd.tsv', columns = ['dim1','dim2'],sep='\t')
        DimRedDF.to_csv('temp/cl.tsv', columns = ['cl'],sep='\t')
        if noEnd:
            cmdToRun= " ".join(["docker run --rm -v", str(Path.cwd())+"/temp/:/data/temp",
                    "slingshot:base /bin/sh -c \"Rscript data/run_slingshot.R",
                    "--input=/data/temp/rd.tsv --input-type=matrix",
                    "--cluster-labels=/data/temp/cl.tsv",
                    "--start-clus="+startClust+'\"'])

        else:
            cmdToRun= " ".join(["docker run --rm -v", str(Path.cwd())+"/temp/:/data/temp",
                                "slingshot:base /bin/sh -c \"Rscript data/run_slingshot.R",
                                "--input=/data/temp/rd.tsv --input-type=matrix",
                                "--cluster-labels=/data/temp/cl.tsv",
                                "--start-clus="+startClust, "--end-clus="+endClust+'\"'])
        print(cmdToRun)
        os.system(cmdToRun)
        for outPath in outPaths:
            os.system("cp temp/PseudoTime.csv "+outPath+"/PseudoTime.csv")
        os.system("rm -rf temp/")

def main(args):
    opts, args = parseArgs(args)
    genSamples(opts)
    
if __name__ == "__main__":
    main(sys.argv)
