import os
from tqdm import tqdm
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

def genDropouts(opts):
    """
    :param outPrefix: default='', Prefix for output files.
    :type outPrefix:  str
    :param expr: Path to expression data file
    :type expr: str
    :param pseudo: Path to pseudotime file
    :type pseudo: str
    :param refNet: Path to reference network file
    :type refNet: str
    :param nCells: default=100, Number of cells to sample.    
    :type nCells: int
    :param dropout: default=False, Carry out dropout analysis?
    :type dropout: bool
    :param drop-cutoff: default=0.5, Specify quantile cutoff on gene expression    
    :type drop-cutoff: float
    :param drop-prob: default=0.5, Specify the probability of dropping a gene below quantile q     
    :type drop-prob: float
    :param samplenum:, default=None, Sample Number
    :type samplenum: int
    """
    
    if opts['dropout']:
        dropoutCutoffs = opts['drop_cutoff']
    else:
        dropoutCutoffs = 0

    ## Read the ExpressionData.csv file
    expDF = pd.read_csv(opts['expr'], index_col=0)
    ## Read the PseudoTime.csv file
    PTDF = pd.read_csv(opts['pseudo'], index_col=0)
    ## Read the refNetwork.csv file
    refDF = pd.read_csv(opts['refNet'], index_col=0)

    ## Generate output path for the dropout datasets
    path = opts['outPrefix'] + '-' +str(opts['nCells'])
    if dropoutCutoffs != 0:
        path += '-' + str(int(100*dropoutCutoffs))+  '-' + str( str(opts['drop_prob']))
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
                    if cointoss < opts['drop_prob']:
                        DropOutDF.loc[idx,col] = 0.0

    DropOutDF.to_csv(path + '/ExpressionData.csv')

def genSamples(opts):
    """
    TODO: Document
    """
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
        timepoints = np.random.choice(range(1,maxtime), size=opts.nCells)
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
        sampledf.to_csv(outfpath + '/ExpressionData.csv')
        ## Read refNetwork.csv
        refdf = pd.read_csv(opts.input_path + '/refNetwork.csv')
        refdf.to_csv(outfpath + '/refNetwork.csv',index=False)
        if len(numclusters) == 1:
            ptdf = pd.DataFrame(np.array(pts),
                            index=pd.Index(cellids),columns = ['PseudoTime'])
        else:
            subsetcluster = clusterdf.loc[[cid.split('_')[0] for cid in cellids]]
            ptdf = pd.DataFrame(index=pd.Index(cellids),
                                columns = ['PseudoTime' + str(1+i) for i in range(len(numclusters))])
            for (cid, row), pt in zip(ptdf.iterrows(), pts):
                ptdf.loc[cid]['PseudoTime' +\
                              str(int(subsetcluster.loc[cid.split('_')[0]]['cl']) + 1)] = round(pt,5)
        ptdf.to_csv(outfpath + '/PseudoTime.csv',na_rep='NA')                        

def doDimRed(opts):
    """
    Carry out dimensionality reduction
    """
    ExpDF = pd.read_csv(opts['expr'],index_col=0, header = 0)
    ptDF = pd.read_csv(opts['pseudo'],index_col=0, header = 0)
    perplexity = opts['perplexity']
    print(perplexity)
    print("Computing TSNE...")
    DimRedRes = TSNE(n_components = 2, perplexity = perplexity).fit_transform(ExpDF.T)
    
    DimRedDF = pd.DataFrame(DimRedRes,columns=['dim1','dim2'],
                            index=pd.Index(list(ExpDF.columns)))
    DimRedDF.loc[:,'pt'] = ptDF.min(axis='columns')    
    DimRedDF.to_csv(str(opts['expr'].parent) + '/tsne' + str(perplexity)+'.tsv', sep='\t')
    plt.figure()
    plt.scatter(DimRedDF.dim1, DimRedDF.dim2, c=DimRedDF.pt)
    plt.savefig(str(opts['expr'].parent)+'/tsne' + str(perplexity) + '.png')
    plt.close()
    
def computeSSPT(opts):
    #nClust, outPath, noEnd = False, perplexity = 500):
    '''
    author: Aditya Pratapa
    Compute PseudoTime using 'slingshot'.
    Needs the input GenexCells expression data frame.
    Needs number of clusters to be expected in the 
    dataset. 
    E.g., Linear: k=1 (no PseudoTime inference is done)
    E.g., Bifurcating: k=3 (1 initial and 2 terminal)
    E.g., Trifurcating: k=4 (1 initial and 3 terminal)
    '''
    ExpDF = pd.read_csv(opts['expr'],index_col=0, header = 0)
    ptDF = pd.read_csv(opts['pseudo'],index_col=0, header = 0)
    nClust = opts['nClusters']
    outPath = opts['outPrefix']
    perplexity = opts['perplexity']
    noEnd = opts['noEnd']
    
    if nClust == 1:
        # Return simulation time as PseduoTime
        ptDF.loc[ExpDF.columns].to_csv(outPath+"/PseudoTime.csv", columns =['Time'])
 
    else:
        ### Compute PseudoTime ordering using slingshot
        
        # Step-1: Compute dimensionality reduction. This is done in doDimRed().
        # Currently only does TSNE
        # TODO: Add PCA

        # Step-2: Read TSNE results to a dataframe
        if not Path(settings['expr'].parent / 'tsne.tsv').is_file():
            # First compute tSNE if this hasn't been done
            doDimRed(opts)
        DimRedDF = pd.read_csv(settings['expr'].parent / 'tsne.tsv',sep='\t')
        
        # Step-3: Compute kMeans clustering
        DimRedDF.loc[:,'cl'] = KMeans(n_clusters = nClust).fit(ExpDF.T).labels_
        # Identify starting cluster

        
        # Step-4: Identify cells corresponding to initial and final states
        # Cells in initial states are identified from cluster with smallest
        # mean experimental time
        # Cells in final states are rest of the clusters
        startClust = DimRedDF.groupby('cl').mean()['pt'].idxmin()
        endClust = ','.join([str(ix) for ix in DimRedDF.groupby('cl').mean()['pt'].index if ix != startClust])
        startClust = str(startClust)
        
        # Step-5: Create a temporary directory and write files necessary to run slingshot
        os.makedirs(outPath, exist_ok = True)

        DimRedDF.to_csv(outPath + '/rd.tsv', columns = ['dim1','dim2'],sep='\t')
        DimRedDF.to_csv(outPath + '/cl.tsv', columns = ['cl'],sep='\t')
        if noEnd:
            cmdToRun= " ".join(["docker run --rm -v", str(Path.cwd()) + "/" +outPath +"/:/data/temp",
                    "slingshot:base /bin/sh -c \"Rscript data/run_slingshot.R",
                    "--input=/data/temp/rd.tsv --input-type=matrix",
                    "--cluster-labels=/data/temp/cl.tsv",
                    "--start-clus="+startClust+'\"'])

        else:
            cmdToRun= " ".join(["docker run --rm -v", str(Path.cwd())+ "/" +outPath +"/:/data/temp",
                                "slingshot:base /bin/sh -c \"Rscript data/run_slingshot.R",
                                "--input=/data/temp/rd.tsv --input-type=matrix",
                                "--cluster-labels=/data/temp/cl.tsv",
                                "--start-clus="+startClust, "--end-clus="+endClust+'\"'])
        print(cmdToRun)
        os.system(cmdToRun)
        # Do this only for the first file
        tn = pd.read_csv(outPath+"/rd.tsv",
             header = 0, index_col =None, sep='\t')

        cl = pd.read_csv(outPath+"/cl.tsv",
                         header = 0, index_col =0, sep='\t')

        curveFile = open(outPath+"/curves.csv","r")
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
        detPT = pd.read_csv(outPath+"/SlingshotPT.csv",
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
        detPT = pd.read_csv(outPath+"/PseudoTime.csv",
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

        tn.to_csv(outPath+"/Updated_rd.tsv",
                  sep='\t')
        plt.savefig(outPath+"/SlingshotOutput.png")
    os.system("rm -rf temp/")    
        
