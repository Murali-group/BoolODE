#Silhouette Model example code below retrieved from scitkit-learn, currently editing to fit/run on BoolODE expression data

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import matplotlib.style as style
from optparse import OptionParser
from pandas import DataFrame

def parseArgs(args):
    parser = OptionParser()
    parser.add_option('', '--expressionfile', type='str', help='Path to ExpressionData.csv file')
    parser.add_option('', '--outPrefix', type='str', default='', help='Prefix for output files')
    
    (opts, args) = parser.parse_args(args)
    return opts, args


def main(args):
    opts, args = parseArgs(args)
    expressionfile = opts.expressionfile
    outPrefix = opts.outPrefix
    if expressionfile is None or len(expressionfile) == 0:
        print("Please specify path to ExpressionData.csv file")
        sys.exit
    if len(expressionfile) > 0:
        expfileDF = pd.read_csv(expressionfile, sep='\t', engine='python', index_col=0)
        data = expfileDF.transpose()
        print(data)
    if len(outPrefix) > 0:
        if '/' in outPrefix:
            outDir = '/'.join(outPrefix.split('/')[:-1])
            if not os.path.exists(outDir):
                print(outDir, "does not exist, creating it...")
                os.makedirs(outDir)
                
        range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10]
        silhouette_avg_n_clusters = []
        
        for n_clusters in range_n_clusters:
            clusterer = KMeans(n_clusters=n_clusters, random_state=42)
            cluster_labels = clusterer.fit_predict(data)
            
            silhouette_avg = silhouette_score(data, cluster_labels)
            print("For n_clusters =", n_clusters, "The average silhouette_score is :", silhouette_avg)
            
            silhouette_avg_n_clusters.append(silhouette_avg)
            
        best_avg_silhouette_value = max(silhouette_avg_n_clusters)
        best_num_cluster = silhouette_avg_n_clusters.index(best_avg_silhouette_value)
        
        print("The best average silhouette score is: ", best_avg_silhouette_value)
        
        #df = pd.DataFrame(silhouette_avg_n_clusters)
        #print(df)
        #df.to_csv(outPrefix + 'silhouettescores.csv')
        silhouette_avg_n_cluster.to_csv(outPrefix + 'silhouettescores.csv')
        
if __name__ == "__main__":
    main(sys.argv)
#Read ExpressionData.csv file
#data = pd.read_csv (r'/Users/cbuck016/Desktop/BoolODE-0.1/VSC_5000/ExpressionData.csv', index_col=0)
#X = data.transpose()
#print(X)

#range_n_clusters = [2, 3, 4, 5, 6, 7, 8, 9, 10]
#silhouette_avg_n_clusters = []

#for n_clusters in range_n_clusters:
    # Initialize the clusterer with n_clusters value and a random generator
    # seed of 10 for reproducibility.
    #clusterer = KMeans(n_clusters=n_clusters, random_state=42)
    #cluster_labels = clusterer.fit_predict(data)

    # The silhouette_score gives the average value for all the samples.
    # This gives a perspective into the density and separation of the formed
    # clusters
    #silhouette_avg = silhouette_score(data, cluster_labels)
    #print("For n_clusters =", n_clusters,
          #"The average silhouette_score is :", silhouette_avg)

    #silhouette_avg_n_clusters.append(silhouette_avg)
    
#best_avg_silhouette_value = max(silhouette_avg_n_clusters)
#best_num_cluster = silhouette_avg_n_clusters.index(best_avg_silhouette_value)

#print("The best average silhouette score is: ", best_avg_silhouette_value)
