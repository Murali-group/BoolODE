import pandas as pd
import numpy as np
import os
from tqdm import tqdm

df = pd.read_csv('ExpressionData.csv',index_col=0)
refdf = pd.read_csv('refNetwork.csv')
seed = 0
np.random.seed(seed)

headers = df.columns
experiments = set([h.split('_')[0] for h in headers])
dfmax = df.values.max()
droplist = []
for e in experiments:
    ecount = 0
    cellcount = 0
    for col in  df.columns:
        if e == col.split('_')[0]:
            cellcount += 1
            if df[col].max() < 0.1*dfmax:
                ecount +=1
    if ecount > 0:
        droplist.append(e)

for d in tqdm(droplist):
    for col in df.columns:
        if d == col.split('_')[0]:
            df.drop(col,axis=1,inplace=True)


headers = df.columns
experiments = list(set([h.split('_')[0] for h in headers]))

subsample = [10,25,50]
numrep = 10

for s in tqdm(subsample):
    for i in range(numrep):
        path = 'output/HSC_' + str(s*100) +'_' + str(i)
        if not os.path.exists(path):
            os.makedirs(path)
        expsInd = np.random.choice(len(experiments),replace=False,size=s)
        newDF = pd.DataFrame(index=df.index)
        for col in df.columns:
            for ei in expsInd:
                if experiments[ei] == col.split('_')[0]:
                    newDF[col] = df[col]
        newDF.to_csv(path +'/ExpressionData.csv')

        cell_id = list(newDF.columns)
        time = [float('.'.join(c.split('_')[1:])) for c in cell_id]
        ptime = [(t - min(time))/(max(time) - min(time)) for t in time]
        exps = []
        for ei in expsInd:
            for _ in range(100):
                exps.append(experiments[ei])
        PseudoTimeDict = {'Cell ID':cell_id, 'PseudoTime':ptime,
                          'Time':time,'Experiment':[c.split('_')[0].replace('E','') for c in cell_id]}
        PseudoTimeDF = pd.DataFrame(PseudoTimeDict)
        PseudoTimeDF[['Cell ID','PseudoTime','Time','Experiment']].to_csv(path +'/PseudoTime.csv',index=False)
        refdf.to_csv(path + '/refNetwork.csv',index=False)

