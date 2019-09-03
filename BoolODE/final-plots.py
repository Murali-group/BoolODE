import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.manifold import TSNE

pathPrefix = '/data/aditya-pratapa/projects/2019-05-01-Benchmarking/RNMethods/inputs/'
simulatedSuffix = '_2000_0/'
dynPrefix = 'dyn-'
dynModels = [
    'LI',
    'LL',
    'CY',
    'BF',
    'BFC',
    'TF'
]

settings = {
    'LI':{'tsne':False},
    'LL':{'tsne':False},
    'BF':{'tsne':False},
    'BFC':{'tsne':False},    
    'CY':{'tsne':False},
    'TF':{'tsne':False}
}
simulatedModelPaths = {'mCAD':'mCAD/mCAD_000_0/',
                       'VSC':'VSC/VSC_2000_0/',
                       'GSD':'GSD/GSD_2000_3/',
                       'HSC':'HSC/HSC_2000_2/'}
simulatedModelFname = 'Updated_rd.tsv'
ptime = '_5000_2'
f, axes = plt.subplots(len(dynModels), 2, figsize=(10, 5*len(dynModels)))
plt.rcParams.update({'font.size': 18})
#for mname,ax in zip(dynModels,axes):
for mname in dynModels:
    print(mname)

    exprDF = pd.read_csv(pathPrefix + dynPrefix + mname + '/'+\
                         dynPrefix + mname + ptime + "/ExpressionData.csv", 
                         header = 0, index_col =0)
    print(exprDF.head())
    if settings[mname]['tsne']:
        dimRed = TSNE(n_components = 2,
                      perplexity=500).fit_transform(exprDF.T)
        
        tn = pd.DataFrame(dimRed, columns=['tSNE-1','tSNE-2'],
                          index=pd.Index(list(exprDF.columns)))
        
        headers = exprDF.columns
        PTData = pd.read_csv(pathPrefix + dynPrefix + mname + '/'+\
                             dynPrefix + mname + ptime +"/PseudoTime.csv", 
                             header = 0, index_col =0)
    
        colNames = PTData.columns
        for idx in range(len(colNames)):
            # Select cells belonging to each pseudotime trajectory
            colName = colNames[idx]
            index = PTData[colName].index[PTData[colName].notnull()]
            tn.loc[index,'cl'] = int(idx)
            #tn['cl'] = tn['cl'].astype('category')
            tn.head()
        
        sns.set_palette("viridis")
        #for col in ptDF.columns:
        #tsneDF[traj].loc[:,col] = ptDF.loc[tsneDF[traj].index,col]
    
        print(tn.cl)
        time = [float('.'.join(t.split('_')[1:])) for t in tn.index]
        tn['time'] = time
        tn.to_csv('dyn-'+mname+'.csv')
        
    tn = pd.read_csv('dyn-'+mname+'.csv',index_col=0)
    f, ax = plt.subplots(1, 1, figsize=(5, 5))    
    sns.scatterplot(x='tSNE-1',y='tSNE-2', 
                    data = tn,
                    hue = 'cl',
                    palette = "Set1",
                    ax = ax) 
    ax.legend_.remove()
    ax.set_ylabel('t-SNE 2')
    ax.set_xlabel('t-SNE 1')
    plt.tight_layout()
    plt.savefig('final-dyn-'+mname+'-cluster.png',dpi=300)    
    f, ax = plt.subplots(1, 1, figsize=(5, 5))    
    sns.scatterplot(x='tSNE-1',y='tSNE-2', 
                    data = tn, hue = 'time',
                    palette = "viridis",ax = ax)
    ax.legend_.remove()
    ax.set_ylabel('t-SNE 2')
    ax.set_xlabel('t-SNE 1')
    plt.tight_layout()
    plt.savefig('final-dyn-'+mname+'-exptime.png',dpi=300)
    
