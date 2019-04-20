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

PC = PCA(n_components=2).fit_transform(Cells)
if tsne_flag:
    embed = TSNE(n_components=2,
                 perplexity=p
    ).fit_transform(PC)
    

colors = [float(h.split('_')[1].replace('|','.')) for h in DF.columns]
experiments = set([h.split('_')[0] for h in DF.columns])

if tsne_flag:
    PCDF = pd.DataFrame(embed,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))
else:
    PCDF = pd.DataFrame(PC,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))
    
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
plt.savefig(inFile.split('/')[1].split('.')[0] + '_dimensionality-reduction.png')
