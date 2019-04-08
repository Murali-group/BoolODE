import pandas as pd
import sys
import numpy as np

inFile = sys.argv[1]
outPrefix = sys.argv[2]


inDF = pd.read_csv(inFile, index_col = 0, header = 0, sep = '\t')
geneDict = {}
geneDict['gene_short_name'] = [gene.replace('x_', '') for gene in inDF.index] # required column!!

Cells = inDF.columns
cellDict = {}
cellDict['Experiment'] = []
cellDict['Time'] = []
for cell in Cells:
    cellDict['Experiment'].append(cell.split('_')[0].replace('E',''))
    cellDict['Time'].append(cell.split('_')[1])

cellDF = pd.DataFrame(cellDict, index = inDF.columns)
cellDF.to_csv(outPrefix + 'Cells.txt', sep = '\t')

geneDF = pd.DataFrame(geneDict, index = inDF.index)
geneDF.to_csv(outPrefix + 'Genes.txt', sep = '\t')


ProcessOutputFlag =  sys.argv[3]

inDFT = inDF.T
if ProcessOutputFlag == 'True':
    print(inDFT.mean())
    inDFCounts = inDFT*np.random.randint(low =1000, high = 1001, size = inDFT.shape)/inDFT.mean()
    print(inDFCounts.head())
    inDFCountsInt = inDFCounts.apply(np.int64)
    inDFCountsInt.T.to_csv(outPrefix + 'ExprMat.txt', sep = '\t')




