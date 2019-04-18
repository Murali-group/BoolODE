from sklearn.manifold import TSNE
#from MulticoreTSNE import MulticoreTSNE as TSNE
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
# names = ['ode',
#          'stoch']
names = ['ode',
         'stoch']
#f, ax = plt.subplots(1,11)
f, ax = plt.subplots(1,len(names))
perp = [30,50]
p = 30
variables = ['x_Gata2','x_EgrNab']
v = 'x_Gata2'

#which = 'tsne'
which = 'pca'

# for n in names:
#     path = 'success-1/tsne'+n+'_experiment.txt'
#     DF = pd.read_csv(path,sep='\t',index_col=0)
#     Cells = DF.T.values


#     PC = PCA(n_components=5).fit_transform(Cells)
#     embed = TSNE(n_components=2,
#                  perplexity=p
#     ).fit_transform(PC)
        
#     #for a,p in zip(arow,perp):
#     variables = DF.index
#     for v,a in zip(variables,ax):
#         #colors = [float(h.split('_')[1]) for h in DF.columns]
#         colors = DF.loc[v].values #[h for h in DF.loc['x_Eklf']]
#         a.scatter(embed[:,0],embed[:,1],c=colors)
#         a.set_title(v)
#         #plt.colorbar()


def plot_colourline(x,y,c):
    c = cm.viridis((c-np.min(c))/(np.max(c)-np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i])
    return

# for n,a in zip(names,ax):
#     path = 'emt_' +n+'_experiment.txt'
#     #path = 'success-1/tsne' +n+'_experiment.txt'    
#     DF = pd.read_csv(path,sep='\t',index_col=0)
#     Cells = DF.T.values
#     PC = PCA(n_components=2).fit_transform(Cells)
#     colors = [float(h.split('_')[1].replace('|','.')) for h in DF.columns]
#     experiments = set([h.split('_')[0] for h in DF.columns])
#     PCDF = pd.DataFrame(PC,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))
#     if which == 'tsne':
#         embed = TSNE(n_components=2,
#                      perplexity=p#,n_jobs=16,
#         ).fit_transform(PC)

#         sc = a.scatter(embed[:,0],embed[:,1],c=colors)
#     else:
#         for e in experiments:
#             toplotX = []
#             toplotY = []
#             colors = []
#             for indexname,row in PCDF.iterrows():
#                 if e in str(indexname):
#                     toplotX.append(row['PC1'])
#                     toplotY.append(row['PC2'])
#                     colors.append(float(indexname.split('_')[1].replace('|','.')))
#             a.plot(toplotX[0],toplotY[0],'ro')                    
#             plot_colourline(a,toplotX,toplotY,colors)

#     a.set_title(n)

n = 'stoch'
path = 'mutin_wpro_' +n+'_experiment.txt'
#path = 'success-1/tsne' +n+'_experiment.txt'    
DF = pd.read_csv(path,sep='\t',index_col=0)
Cells = DF.T.values
PC = PCA(n_components=2).fit_transform(Cells)
embed = TSNE(n_components=2,
             perplexity=p#,n_jobs=16,
).fit_transform(PC)





colors = [float(h.split('_')[1].replace('|','.')) for h in DF.columns]
experiments = set([h.split('_')[0] for h in DF.columns])
#PCDF = pd.DataFrame(PC,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))
PCDF = pd.DataFrame(embed,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))

for e in experiments:
    toplotX = []
    toplotY = []
    colors = []
    for indexname,row in PCDF.iterrows():
        if e in str(indexname):
            toplotX.append(row['PC1'])
            toplotY.append(row['PC2'])
            colors.append(float(indexname.split('_')[1].replace('|','.')))
    plt.plot(toplotX[0],toplotY[0],'ro')
    plt.scatter(toplotX,toplotY,c=colors)                    
    #plot_colourline(toplotX[:-1],toplotY[:-1],colors)

    
plt.legend()
plt.show()
