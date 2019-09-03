import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from optparse import OptionParser
import pandas as pd
import seaborn as sns
import os
parser = OptionParser()

parser.add_option('-m','--model',default='',type=str,
                  help='Specify model name')

(opts, args) = parser.parse_args()
modelName = opts.model
if not os.path.exists(modelName +'/'):
    os.mkdir(modelName)


tsnedf = pd.read_csv(modelName + '_ExpressionData-dropped.csv_dimred.txt',sep=',',index_col=0)
expdf = pd.read_csv(modelName + '_ExpressionData-dropped.csv',sep=',',index_col=0)

sizeDict = {'mCAD':(2,3),
            'VSC':(3,3),
            'HSC':(3,4),
            'GSD_input':(4,5)}

axrows = sizeDict[modelName][0]
axcols = sizeDict[modelName][1]
length = 4*axcols

#f = plt.figure(figsize=(length, length*float(axrows)/float(axcols)))

mycmap = sns.cubehelix_palette(as_cmap=True)
oog = {'mCAD':['Sp8','Fgf8','Pax6',
                'Coup','Emx2'],
       'VSC':['Irx3','Pax6','Nkx62',
              'Nkx22','Nkx61','Olig2',
              'Dbx2','Dbx1']}
if modelName in oog:
    orderofgenes = oog[modelName]
else:
    orderofgenes = list(expdf.index)
#for axcount,(gene,row) in enumerate(expdf.iterrows()):
for axcount,gene in enumerate(orderofgenes):
    f = plt.figure(figsize=(4.5,4))
    ax = f.gca()
    #ax = f.add_subplot(axrows, axcols, axcount+1)
    geneexpval = expdf.loc[gene].values
    sns.scatterplot(x='tsne1',y='tsne2',data=tsnedf,hue=geneexpval,ax=ax)
    ax.set_title(gene)
    norm = plt.Normalize(min(geneexpval), max(geneexpval))
    sm = plt.cm.ScalarMappable(norm=norm,cmap=mycmap)
    sm.set_array([])
    ax.get_legend().remove()
    ax.figure.colorbar(sm)
    ax.set_xlabel('tSNE1')
    ax.set_ylabel('tSNE2')    
    plt.tight_layout()
    plt.savefig(modelName + '/' + gene + '.png',dpi=200) 
#f.suptitle(modelName)
# plt.tight_layout()
# plt.savefig(modelName + '_gene.png',dpi=200) 
#plt.savefig(modelName + '_gene.pdf')
    
