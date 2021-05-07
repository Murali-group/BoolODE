# Analysis of silhouette scores and best predicted number of clusters for 1,000 BoolODE simulations 

import csv
my_reader = csv.reader(open('/home/cbuck016/BoolODE-0.1/mCAD-sims/mCAD-silhouettescores.csv'))
ctr = 0
for cluster in my_reader:
  if cluster[1] == '2':
    ctr += 1
print(ctr)
                       
