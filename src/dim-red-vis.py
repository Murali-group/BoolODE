from sklearn.manifold import TSNE
#from MulticoreTSNE import MulticoreTSNE as TSNE
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from optparse import OptionParser 
parser = OptionParser()

parser.add_option('-i','--inFile',default='',type=str,
                  help='Specify input expression matrix file name')
parser.add_option('-t','--tsne',action='store_true',default=False,
                  help='Visualized tsne instead of default PCA')

(opts, args) = parser.parse_args()
inFile = opts.inFile
tsne_flag = opts.tsne

f, ax = plt.subplots(1,1)

which = 'pca'


def plot_colourline(x,y,c):
    c = cm.viridis((c-np.min(c))/(np.max(c)-np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=c[i])
    return


n = 'stoch'

DF = pd.read_csv(inFile,sep=',',index_col=0)
Cells = DF.T.values

####################
# Do PCA and tSNE
PC = PCA(n_components=2).fit_transform(Cells)
embed = TSNE(n_components=2).fit_transform(Cells)
####################    

colors = [float(h.split('_')[1]) for h in DF.columns]
experiments = set([h.split('_')[0] for h in DF.columns])
PCDF = pd.DataFrame(PC,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))

PCDF['tsne1'] = embed[:,0]
PCDF['tsne2'] = embed[:,1]
    
PCDF['time'] = colors
PCDF.to_csv(inFile + '_dimred.txt')

f,ax = plt.subplots(2,1)
for e in experiments:
    toplotX = []
    toplotY = []
    colors = []
    for indexname,row in PCDF.iterrows():
        if e in str(indexname):
            toplotX.append(row['PC1'])
            toplotY.append(row['PC2'])
            colors.append(float(indexname.split('_')[1]))#.replace('|','.')))
    #ax[0].plot(toplotX[0],toplotY[0],'ro')
    ax[0].scatter(toplotX,toplotY,c=colors)
    ax[0].set_title('PCA')
    for indexname,row in PCDF.iterrows():
        if e in str(indexname):
            toplotX.append(row['tsne1'])
            toplotY.append(row['tsne2'])
            colors.append(float(indexname.split('_')[1]))#.replace('|','.')))
    #ax[1].plot(toplotX[0],toplotY[0],'ro')
    ax[1].scatter(toplotX,toplotY,c=colors)
    ax[1].set_title('tSNE')
    #plot_colourline(toplotX[:-1],toplotY[:-1],colors)
    
plt.legend()
plt.savefig(inFile.split('.csv')[0] + '_dimensionality-reduction.png')
