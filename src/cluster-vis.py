import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-i','--inFile',default='',type=str,
                  help='Specify input expression matrix file name')
(opts, args) = parser.parse_args()
inFile = opts.inFile
DF = pd.read_csv(inFile,sep=',',index_col=0)
genes = DF.loc[[row for row in DF.index if 'p_' not in row]]
print(genes.shape)
genes.drop([col for col in genes.columns if genes[col].sum() < 0.2*len(genes.index)],axis=1,inplace=True)
sns.clustermap(genes)
plt.savefig(inFile.split('.')[0],'_clustered-genes.png')
