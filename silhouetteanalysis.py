# Analysis of silhouette scores and best predicted number of clusters for 1,000 BoolODE simulations 

import csv
my_reader = csv.reader(open('/home/cbuck016/BoolODE-0.1/mCAD-sims/mCAD-silhouettescores.csv'))
ctr_of_2 = 0
ctr_of_3 = 0
ctr_of_4 = 0
ctr_of_5 = 0
ctr_of_6 = 0
ctr_of_7 = 0
ctr_of_8 = 0
ctr_of_9 = 0
ctr_of_10 = 0

#count the number of times the steady state is reached in the model
for cluster in my_reader:
  if cluster[1] == '2':
    ctr_of_2 += 1
                        
  for cluster in my_reader:
    if cluster[1] == '3':
      ctr_of_3 += 1
      
    for cluster in my_reader:
      if cluster[1] == '4':
        ctr_of_4 += 1
                        
      for cluster in my_reader:
        if cluster[1] == '5':
          ctr_of_5 += 1
          
        for cluster in my_reader:
          if cluster[1] == '6':
            ctr_of_6 += 1
            
          for cluster in my_reader:
            if cluster[1] == '7':
              ctr_of_7 += 1
                   
            for cluster in my_reader:
              if cluster[1] == '8':
                ctr_of_8 += 1
                
              for cluster in my_reader:
                if cluster[1] == '9':
                  ctr_of_9 += 1
                  
                for cluster in my_reader:
                  if cluster[1] == '10':
                    ctr_of_10 += 1
                    
                 
print("2 clusters predicted " + str(ctr_of_2) + " times") 
print("3 clusters predicted " + str(ctr_of_3) + " times")
print("4 clusters predicted " + str(ctr_of_4) + " times") 
print("5 clusters predicted " + str(ctr_of_5) + " times")                 
print("6 clusters predicted " + str(ctr_of_6) + " times")                
print("7 clusters predicted " + str(ctr_of_7) + " times")
print("8 clusters predicted " + str(ctr_of_8) + " times")
print("9 clusters predicted " + str(ctr_of_9) + " times")
print("10 clusters predicted " + str(ctr_of_10) + " times")
