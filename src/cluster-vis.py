import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-p','--prefix',default='',type=str,
                  help='Specify file name prefix [prefix]stoch_experiment.txt')
(opts, args) = parser.parse_args()
prefix = opts.prefix
DF = pd.read_csv(prefix + 'stoch_experiment.txt',sep='\t',index_col=0)
genes = DF.loc[[row for row in DF.index if 'x_' in row]]
genes.drop([col for col in genes.columns if genes[col].sum() < 0.2*len(genes.index)],axis=1,inplace=True)
sns.clustermap(genes)
plt.show()
