# This is a modified version of genVis.py authored by Jon Mallen. Original genVis.py script authored by Murali et al.
# found at https://github.com/Murali-group/BoolODE.

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import itertools

parser = argparse.ArgumentParser("Visualize the BoolODE simulation of cells")
parser.add_argument('-p', '--pathToFiles', default='', type=str, help='Specify path to files')
parser.add_argument('-t', '--tsne', action='store_true', default=False, help='Visualize tsne instead of default PCA')
parser.add_argument('-d', '--dimensions', default=2, type=int, help='Choose dimensions for t-SNE (2 or 3)')
parser.add_argument('-c', '--handleExternalClusterFile', action='store_true', default=False,
                    help='Call to handle an external cluster file')
parser.add_argument('-n', '--plotTitle', default='', nargs='*', help='Write a plot title')

args = parser.parse_args()
path = args.pathToFiles
inFile = path + "/expressiondata.csv"
timeFile = path + "/pseudotime.csv"
tsne_flag = args.tsne
tsne_dim = args.dimensions
clusterFile = args.pathToFiles + "/clusterids.csv"

DF = pd.read_csv(inFile, sep=',', index_col=0)
Cells = DF.T.values

####################
# Do PCA and tSNE
PC = PCA(n_components=2).fit_transform(Cells)
embed = TSNE(n_components=tsne_dim).fit_transform(Cells)
####################
ptDF = pd.read_csv(timeFile, sep=',', index_col=0)
color_scale = max(ptDF["Time"])
colors_raw = [h.split('_')[1] for h in DF.columns]
colors_raw = [int(i) for i in colors_raw]
colors = [x / color_scale for x in colors_raw]
experiments = set([h.split('_')[0] for h in DF.columns])
PCDF = pd.DataFrame(PC, columns=['PC1', 'PC2'], index=pd.Index(list(DF.columns)))
PCDF['tsne1'] = embed[:, 0]
PCDF['tsne2'] = embed[:, 1]
if tsne_dim == 3:
    PCDF['tsne3'] = embed[:, 2]

PCDF['time'] = colors
PCDF.to_csv(inFile + '_dimred.txt')

# TODO: implement graphing option for PCA

f, ax = plt.subplots(1, 2, figsize=(10, 5))
if tsne_dim == 3:
    ax[0].set_axis_off()
    ax[0] = f.add_subplot(1, 2, 1, projection="3d")
    ax[0].scatter3D(PCDF['tsne1'], PCDF['tsne2'], PCDF['tsne3'], c=PCDF['time'])
    ax[0].set_zlabel('t-SNE 3')
else:
    ax[0].scatter(PCDF['tsne1'], PCDF['tsne2'], c=PCDF['time'])
ax[0].set_xlabel('t-SNE 1')
ax[0].set_ylabel('t-SNE 2')
ax[0].set_aspect('auto')
ax[0].set_title('Simulation Time')

if args.handleExternalClusterFile:
    CF = pd.read_csv(clusterFile, sep=',', index_col=0)
    cluster_colors_raw = CF['cl'].tolist()
    cluster_color_scale = max(CF['cl'])
    cluster_colors = [y / color_scale for y in cluster_colors_raw]
else:
    cluster_colors = list(itertools.repeat(.5, len(DF.columns)))

plt.rcParams['image.cmap'] = 'Spectral'
if tsne_dim == 3:
    ax[1].set_axis_off()
    ax[1] = f.add_subplot(1, 2, 2, projection="3d")
    ax[1].scatter3D(PCDF['tsne1'], PCDF['tsne2'], PCDF['tsne3'], c=cluster_colors)
    ax[1].set_zlabel('t-SNE 3')
else:
    ax[1].scatter(PCDF['tsne1'], PCDF['tsne2'], c=cluster_colors)
ax[1].set_xlabel('t-SNE 1')
ax[1].set_ylabel('t-SNE 2')
ax[1].set_aspect('auto')
ax[1].set_title('Clusters')

plot_title = ' '.join(args.plotTitle)
plt.suptitle(plot_title, fontsize=20)

if tsne_dim == 3:
    plt.show()
else:
    plt.savefig(inFile.split('.csv')[0] + '_dimensionality-reduction.png')
