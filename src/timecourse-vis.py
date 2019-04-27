import matplotlib.pyplot as plt
import pandas as pd
import sys
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-i','--inFile',default='',type=str,
                  help='Specify input expression matrix file name')
(opts, args) = parser.parse_args()
inFile = opts.inFile

DF = pd.read_csv(inFile, sep=',',index_col=0)
genes = DF.index
f,ax = plt.subplots(len(genes),1,figsize=(2*len(genes),10))

tp = DF.columns
experiments = set([t.split('_')[0] for t in tp])
times = sorted(list(set([float(t.split('_')[1]+'.'+t.split('_')[2]) for t in tp])))
for g,a in zip(genes,ax):
    V = DF.loc[g].values
    for e in experiments:
        toplot = []
        for h,v in zip(tp,V):
            if e+'_' in h:
                toplot.append(v)
        a.plot(times,toplot)
    a.set_ylabel(g.replace('x_',''))        

a.set_xlabel('STOCH')
plt.savefig(inFile.split('.csv')[0] + '_time-course.png')
