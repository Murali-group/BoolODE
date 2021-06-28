# This is a modified version of genVis.py authored by Jon Mallen. Original genVis.py script authored by Murali et al.
# found at https://github.com/Murali-group/BoolODE.

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import itertools

parser = argparse.ArgumentParser("Visualize the BoolODE simulation of cells")
parser.add_argument('-f', '--pathToFiles', default='', type=str, help='Specify path to files')
parser.add_argument('-p', '--pca', nargs='*', help='Call to visualize with PCA, specify dimension (2 or 3) as argument')
parser.add_argument('-t', '--tsne', nargs='*', help='Call to visualize with t-SNE, '
                                                    'specify dimension (2 or 3) as argument')
parser.add_argument('-u', '--umap', nargs='*', help='Call to visualize with UMAP, '
                                                    'specify dimension (2 or 3) as argument')
parser.add_argument('-c', '--handleExternalClusterFile', action='store_true', default=False,
                    help='Call to handle an external cluster file')
parser.add_argument('-n', '--plotTitle', default='', nargs='*', help='Write a plot title')

# Parse arguments
args = parser.parse_args()
path = args.pathToFiles
inFile = path + "/expressiondata.csv"
timeFile = path + "/pseudotime.csv"
pca_flag = args.pca is not None
tsne_flag = args.tsne is not None
umap_flag = args.umap is not None
clusterFile = args.pathToFiles + "/clusterids.csv"

# Do PCA, tSNE, UMAP
DF = pd.read_csv(inFile, sep=',', index_col=0)
Cells = DF.T.values
PCDF = pd.DataFrame(index=pd.Index(list(DF.columns)))
if pca_flag:
    if len(args.pca) == 0:
        pca_dim = 2
    else:
        pca_dim = int(args.pca[0])
    PC = PCA(n_components=pca_dim).fit_transform(Cells)
    PCDF['PC1'] = PC[:, 0]
    PCDF['PC2'] = PC[:, 1]
    if pca_dim == 3:
        PCDF['PC3'] = PC[:, 2]
if tsne_flag:
    if len(args.tsne) == 0:
        tsne_dim = 2
    else:
        tsne_dim = int(args.tsne[0])
    tsne_embed = TSNE(n_components=tsne_dim).fit_transform(Cells)
    PCDF['tsne1'] = tsne_embed[:, 0]
    PCDF['tsne2'] = tsne_embed[:, 1]
    if tsne_dim == 3:
        PCDF['tsne3'] = tsne_embed[:, 2]
if umap_flag:
    if len(args.umap) == 0:
        umap_dim = 2
    else:
        umap_dim = int(args.umap[0])
    umap_embed = umap.UMAP(n_components=umap_dim).fit_transform(Cells)
    PCDF['umap1'] = umap_embed[:, 0]
    PCDF['umap2'] = umap_embed[:, 1]
    if umap_dim == 3:
        PCDF['umap3'] = umap_embed[:, 2]


# Prepare time-dependent color scheme
ptDF = pd.read_csv(timeFile, sep=',', index_col=0)
color_scale = max(ptDF["Time"])
colors_raw = [h.split('_')[1] for h in DF.columns]
colors_raw = [int(i) for i in colors_raw]
colors = [x / color_scale for x in colors_raw]
experiments = set([h.split('_')[0] for h in DF.columns])
PCDF['time'] = colors

# Prepare cluster-dependent color scheme
if args.handleExternalClusterFile:
    CF = pd.read_csv(clusterFile, sep=',', index_col=0)
    cluster_colors_raw = CF['cl'].tolist()
    cluster_color_scale = max(CF['cl'])
    cluster_colors = [y / color_scale for y in cluster_colors_raw]
else:
    cluster_colors = list(itertools.repeat(.5, len(DF.columns)))
# PCDF = pd.DataFrame(PC, columns=['PC1', 'PC2'], index=pd.Index(list(DF.columns)))
# PCDF['tsne1'] = embed[:, 0]
# PCDF['tsne2'] = embed[:, 1]
# if tsne_dim == 3:
#     PCDF['tsne3'] = embed[:, 2]

# Write data to text file
PCDF.to_csv(inFile + '_dimred.txt')

# Graphing
plot_title = ' '.join(args.plotTitle)

# PCA graphing


# Common method for plotting

def make_subplot(axis1, axis1label, axis2, axis2label, axis3, axis3label, dim):
    f, ax = plt.subplots(1, 2, figsize=(10, 5))
    plt.rcParams['image.cmap'] = 'viridis'
    if dim == 3:
        ax[0].set_axis_off()
        ax[0] = f.add_subplot(1, 2, 1, projection="3d")
        ax[0].scatter3D(PCDF[axis1], PCDF[axis2], PCDF[axis3], c=PCDF['time'])
        ax[0].set_zlabel(axis3label)
    else:
        ax[0].scatter(PCDF[axis1], PCDF[axis2], c=PCDF['time'])
    ax[0].set_xlabel(axis1label)
    ax[0].set_ylabel(axis2label)
    ax[0].set_aspect('auto')
    ax[0].set_title('Simulation Time')

    plt.rcParams['image.cmap'] = 'Spectral'
    if dim == 3:
        ax[1].set_axis_off()
        ax[1] = f.add_subplot(1, 2, 2, projection="3d")
        ax[1].scatter3D(PCDF[axis1], PCDF[axis2], PCDF[axis3], c=cluster_colors)
        ax[1].set_zlabel(axis3label)
    else:
        ax[1].scatter(PCDF[axis1], PCDF[axis2], c=cluster_colors)
    ax[1].set_xlabel(axis1label)
    ax[1].set_ylabel(axis2label)
    ax[1].set_aspect('auto')
    ax[1].set_title('Clusters')

    plt.suptitle(plot_title, fontsize=20)


# t-SNE graphing
if tsne_flag:
    make_subplot('tsne1', 't-SNE 1', 'tsne2', 't-SNE 2', 'tsne3', 't-SNE 3', tsne_dim)
    plt.savefig(inFile.split('.csv')[0] + '_tSNE_%sd.png' % tsne_dim)

# PCA graphing
if pca_flag:
    make_subplot('PC1', 'PCA 1', 'PC2', 'PCA 2', 'PC3', 'PCA 3', pca_dim)
    plt.savefig(inFile.split('.csv')[0] + '_PCA_%sd.png' % pca_dim)

# UMAP graphing
if umap_flag:
    make_subplot('umap1', 'UMAP 1', 'umap2', 'UMAP 2', 'umap3', 'UMAP 3', umap_dim)
    plt.savefig(inFile.split('.csv')[0] + '_UMAP_%sd.png' % umap_dim)

plt.show()
