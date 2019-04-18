import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
DF = pd.read_csv('test_stoch_experiment.txt',sep='\t',index_col=0)
genes = DF.loc[[row for row in DF.index if 'x_' in row]]
#genes.drop([col for col in genes.columns if genes[col].sum() < 0.2*len(genes.index)],axis=1,inplace=True)
sns.clustermap(genes)
plt.show()
