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

parser.add_option('-p','--pseudoTimeFile',default='',type=str,
                  help='Specify pseudoTimeFile file name')

parser.add_option('-t','--tsne',action='store_true',default=False,
                  help='Visualized tsne instead of default PCA')

(opts, args) = parser.parse_args()
inFile = opts.inFile
tsne_flag = opts.tsne




DF = pd.read_csv(inFile,sep=',',index_col=0)
Cells = DF.T.values

####################
# Do PCA and tSNE
PC = PCA(n_components=2).fit_transform(Cells)
embed = TSNE(n_components=2).fit_transform(Cells)
####################    
ptDF = pd.read_csv(opts.pseudoTimeFile, sep=',', index_col=0)

colors = ptDF.min(axis='columns').values
print(colors)
experiments = set([h.split('_')[0] for h in DF.columns])
PCDF = pd.DataFrame(PC,columns=['PC1','PC2'],index=pd.Index(list(DF.columns)))

PCDF['tsne1'] = embed[:,0]
PCDF['tsne2'] = embed[:,1]
    
PCDF['time'] = colors
PCDF.to_csv(inFile + '_dimred.txt')

f,ax = plt.subplots(2,1,figsize=(5,10))

colors = [h.split('_')[1] for h in DF.columns]
ax[0].scatter(PCDF['PC1'], PCDF['PC2'], c= PCDF['time'])
ax[0].set_title('PCA')
ax[1].scatter(PCDF['tsne1'], PCDF['tsne2'], c= PCDF['time'])
ax[1].set_title('tSNE')

    
plt.legend('')
plt.savefig(inFile.split('.csv')[0] + '_dimensionality-reduction.png')    
