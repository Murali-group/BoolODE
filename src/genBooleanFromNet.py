import pandas as pd
import sys
def main():
    inFile = sys.argv[1]
    outFile = sys.argv[2]
    
    print("Reading network from ", inFile)
    inDF = pd.read_csv(inFile, sep = '\t')
    
    outDict = {}
    outDict['Gene'] = inDF['Gene1'].unique()
    outDict['Rule'] = []
    for gene in outDict['Gene']:
        act = []
        rep = []
        for idx,row in inDF.loc[inDF['Gene1'] == gene].iterrows():
            if row['Type'] == '+':
                act.append(row['Gene2'])
            else:
                rep.append(row['Gene2'])
        if len(act) != 0 and len(rep) != 0:
        	finalRule = '(( ' + ' or '.join(act) + ' ) and not( ' + ' or '.join(rep) + ' ))'
        elif len(act) == 0:
            finalRule = 'not ( ' + ' or '.join(rep) + ' )'
        else:
            finalRule = '( ' + ' or '.join(act) + ' )'
        outDict['Rule'].append(finalRule)
    
    outDF = pd.DataFrame(outDict)
    print(outDF.head())
    print("Writing file to ", outFile)
    outDF.to_csv(outFile, sep = '\t', index = False)

if __name__ == '__main__':
    main()
    

