#Silhouette Method for use with output results from BoolODE datasets -- IN PROGRESS

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import matplotlib.style as style

pathtoBoolFiles = '/home/cbuck016/BoolODE-0.1/mCAD-sims/'
model_name = 'mCAD'
experiment_name = model_name + '-ts-' + str(ts)
pathtoExpressionDataFile = pathtoBoolFiles + '

filenames = ['sim-1-ExpressionData.csv', 'sim-2-ExpressionData.csv'...]
dataframes = []
for f in filenames:
  dataframes.append(pd.read_csv(f))
