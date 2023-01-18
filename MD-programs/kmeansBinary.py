#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:43:52 2022

@author: abnerkahansmack
"""
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import NMF
from scipy.spatial.distance import cdist
import sklearn.cluster
from sklearn.model_selection import train_test_split
from sklearn.metrics import silhouette_score


def positionmap(file, residues):
    distancefile = open(file)
    reader =  distancefile.read()
    
    superlist = reader.split ('} {')
    #print(superlist [-50:])
    #print("\n\n\n")
    superlist [0] = superlist[0].replace('{', '')
    superlist [-1] = superlist[-1].replace('}', '')
    superlist [-1] = superlist[-1][:-1]
    #print(superlist [-50:])
    #print(len(superlist))
    #print(superlist[-9:])
    distancefile.close()
     
    frames = int(len(superlist)/residues)
    contactNP = np.zeros((residues ** 2, frames ))
    reset  =  residues ** 2                  
    counter = 0
    fra = 0                   
    for entry in superlist:
        entry = entry.replace(" ", "")
        #print(entry)
        if  counter < reset:
            for bob in entry:
                #print(bob, counter)
                contactNP [counter,fra] = int(bob)
                counter +=1
        else:
            counter = 0
            fra +=1
            for bob in entry:
                #print(bob, counter)
                contactNP [counter,fra] = int(bob)
                counter +=1


    
    #customW = np.ones((residues**2,numNMF))      
   
    zerolist =[]
    for i in range(residues):
        zerolist.append(i*residues+i)
        
    G2bonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8)]
    G2bondsB = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8), (1,6), (2,4), (1,3), (2,7)]
    G2fbonds = [(0,1), (0,2), (2,3), (3,4), (4,5), (5,6), (3,7),  (7,8), (8,9)  ]
    G2fbondsB = [(0,1), (0,2), (2,3), (3,4), (4,5), (5,6), (3,7),  (7,8), (8,9), (1,2), (2,4), (3,8), (3,5), (2,7)  ]
    M9bonds = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5)]
    M9bondsC = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5), (1,3), (3,5),  (8,10), (2,6), (1,8), (3,7), (2,9), (2,10) ]
    M9bondsB = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5), (1,3),(3,5), (8, 10), (2,6), (1,8), (3,7), (2,9), (2,4) ]
    N2bonds = [(0,1), (1,2), (2,3), (3,4), (2,5), (5,6)]
    N2bondsB = [(0,1), (1,2), (2,3), (3,4), (2,5), (5,6),(1,5), (2,4), (2,6), (1,3)]
    N2fbonds = [(0,1), (0,2), (2,3), (3,6), (6,7), (3,4), (4,5)   ]
    N2fbondsB = [(0,1), (0,2), (2,3), (3,6), (6,7), (3,4), (4,5), (1,2), (3,7), (2,6) , (2,4), (3,5)  ]
    S2bonds = [(0,1), (1,2), (2,3), (3,4),(4,5), (5,6), (2,7), (7,8), (8,9), (9,10)  ]
    S2bondsB = [(0,1), (1,2), (2,3), (3,4),(4,5), (5,6), (2,7), (7,8), (8,9), (9,10), (1,7), (8,10), (4,6), (2,4), (2,8), (1,3)  ]
    
    for bond in G2bondsB:
        zerolist.append((bond[0]*residues)+bond[1])
        zerolist.append((bond[1]*residues)+bond[0])
    print(zerolist)
    for frame in range(frames):
        for zero in zerolist:
            contactNP[zero,frame] = 0
            #if frame ==0:
            #    customW[zero,0] =0
              #  if numNMF == 2 :
                #    customW[zero,1] =0
    return contactNP
def calculate_WSS(points, kmax):
  sse = []
  info = []
  silhouette = []
  print(len(points))
  for k in range(2, kmax+1):
    kmeans = sklearn.cluster.KMeans(n_clusters = k, tol=1**-12).fit(points)
    #centroids = kmeans.cluster_centers_
    #print(centroids, "\n")
    pred_clusters = kmeans.predict(points)
    #print (pred_clusters)
    info.append(pred_clusters)
    
    sse.append(np.sum( kmeans.inertia_))
    silhouette.append(silhouette_score(points, pred_clusters))
    
    
    # calculate square of Euclidean distance of each point from its cluster center and add to current WSS
    #for i in range(len(points)):
      #curr_center = centroids[pred_clusters[i]]
      #curr_sse += (points[i, 0] - curr_center[0]) ** 2 + (points[i, 1] - curr_center[1]) ** 2
      
  return sse, info, silhouette
pm =positionmap('vmd/g2-45d.txt',9 ).transpose()
trainData, testData = train_test_split(pm)
yellow = calculate_WSS(pm, 20)
def findVals(clusters,index, range1, range2, value):
    info = clusters[1][index][range1:range2]
    info = np.where(info ==value)
    print(info)
    
    
print(yellow[0])      
plt.title("Means squared error")
plt.plot(range(2,21),yellow[0], 'bx-')
plt.xticks(range(2,21,2))
plt.xlabel("Number of clusters") 
plt.show()
plt.clf()

plt.title("Silhouette squared error")
plt.plot(range(2,21),yellow[2], 'rx-')
plt.xticks(range(2,21,2))
plt.xlabel("Number of clusters") 
plt.show()
plt.clf()
            
#np.random.choice(np.where(yellow[1][6] == 0)[0])
#don't flatten out matrix