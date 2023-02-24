#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:17:56 2022

@author: abnerkahansmack
"""


import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import NMF
from scipy.spatial.distance import cdist
import sklearn.cluster
from scipy.signal import find_peaks

def positionmap(file, residues):
    mapGL2 =  open(file, 'r')
    maptext = mapGL2.read()  
#print(maptext)
    mapGL2.close()
#print(maptext)
#print(type(maptext))
#print(type(maptext[0]))
#print(maptext.count('\n'))
    contactList = maptext.split("} {")
    contactList [0] = contactList[0].replace('{', '')
    contactList [-1] = contactList[-1].replace('}', '')
   # print(contactList[:50])
    print(len(contactList))
    Nframes = len(contactList)/ residues
    
    
    contactNP = np.zeros((3, residues,int(len(contactList)/residues) ) )
    position = 0
    frame = 0
    for stringy in contactList:
        if position  == residues:
            frame += 1
            position = 0
        contactNP[:, position, frame] = stringy.split(" ")
        position += 1
    #print (contactNP[:,:, :3])
    
    return contactNP 
    sklearn.cluster.SpectralClustering(n_clusters = 8).fit(contactNP)
def moveAverage(L,n):
    return np.convolve(L, np.ones(n), 'valid') / n
def distEu(triple1,triple2):
    return np.sqrt(( (triple1[0] - triple2[0])**2 + (triple1[1] - triple2[1])**2 + (triple1[2] - triple2[2])**2 ))
def distancemap(positionmap):
    dims = positionmap.shape
    distanceNP = np.zeros((dims[1]**2, dims[2]))
    for frame in range(dims[2]):
       position = 0   
       for resid2 in range (dims[1]):
           for resid1 in range (dims[1]):
               
                distanceNP[position,frame] = distEu(positionmap[:, resid2, frame], positionmap[:, resid1, frame])
                position += 1
    return distanceNP.transpose()
        #print (numString)
def calculate_WSS(points, kmax):
  sse = []
  info = []
  print(len(points))
  for k in range(1, kmax+1):
    kmeans = sklearn.cluster.KMeans(n_clusters = k, tol=1**-12).fit(points)
    #centroids = kmeans.cluster_centers_
    #print(centroids, "\n")
    pred_clusters = kmeans.predict(points)
    
    info.append(pred_clusters)
    
    sse.append(np.sum( kmeans.inertia_))
    
    
    # calculate square of Euclidean distance of each point from its cluster center and add to current WSS
    #for i in range(len(points)):
      #curr_center = centroids[pred_clusters[i]]
      #curr_sse += (points[i, 0] - curr_center[0]) ** 2 + (points[i, 1] - curr_center[1]) ** 2
      
  return sse, info

mappo = distancemap(positionmap('vmd/g2-center.txt', 9))
wss = calculate_WSS(mappo, 10)
print(wss)
#tom2 = moveAverage(tom, 10)
#peaks = find_peaks(tom2*-1, rel_height=.3)
#print(peaks)

