import matplotlib.pyplot as plt
import pandas as pd
import sys

#pre = 'success-1/tsne'
pre = ''

odefn = 'ode_experiment.txt'
stochfn = 'stoch_experiment.txt'

#DF = pd.read_csv(pre + odefn, sep='\t',index_col=0)
DF = pd.read_csv(pre + stochfn, sep='\t',index_col=0)
genes = DF.index
genes = genes#[0:-1]
f,ax = plt.subplots(len(genes),2,figsize=(10,10))
tp = DF.columns
experiments = set([t.split('_')[0] for t in tp])
print(experiments)
# for g,a in zip(genes,ax):
#     V = DF.loc[g].values
#     for e in experiments:
#         toplot = []
#         for h,v in zip(tp,V):
#             if e+'_' in h:
#                 toplot.append(v)
#         a[0].plot(toplot)

#     a[0].set_ylabel(g.replace('x_',''))



#DF = pd.read_csv(pre + stochfn, sep='\t',index_col=0)
for g,a in zip(genes,ax):
    V = DF.loc[g].values

    for e in experiments:
        toplot = []
        for h,v in zip(tp,V):
            if e+'_' in h:
                toplot.append(v)
        print(len(toplot))
        a[1].plot(toplot)

    a[1].set_ylabel(g.replace('x_',''))        

a[0].set_xlabel('ODE')
a[1].set_xlabel('STOCH')
#plt.show()
plt.savefig('test.png')
