#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 14:43:59 2022

@author: abnerkahansmack
"""

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import NMF
from scipy.spatial.distance import cdist
#cdist

def distEu(triple1,triple2):
    return np.sqrt(( (triple1[0] - triple2[0])**2 + (triple1[1] - triple2[1])**2 + (triple1[2] - triple2[2])**2 ))
def positionmap(file, residues):
    distancefile = open(file)
    reader =  distancefile.read()
    superlist = reader.split ('} {')
    superlist [0] = superlist[0].replace('{', '')
    superlist [-1] = superlist[-1].replace('}', '')
    distancefile.close()
    contactNP = np.zeros((3, residues,int(len(superlist)/residues) ) )
    position = 0
    frame = 0
    for stringy in superlist:
        if position  == residues:
            frame += 1
            position = 0
        contactNP[:, position, frame] = stringy.split(" ")
        position += 1
    return contactNP 

#g2 = positionmap('vmd/g2-center.txt', 9)        
#for  n in range(9): 
 #   plt.plot(positionmap('vmd/g2-center.txt', 9)[2, n, ::100])
def moveAverage(L,n):
    return np.convolve(L, np.ones(n), 'valid') / n
def distancemap(positionmap):
    dims = positionmap.shape
    print(dims)
    distanceNP = np.zeros((dims[1]**2, dims[2]))
    
    for frame in range(dims[2]):
       position = 0   
       for resid2 in range (dims[1]):
           for resid1 in range (dims[1]):
               
                distanceNP[position,frame] = distEu(positionmap[:, resid2, frame], positionmap[:, resid1, frame])
                position += 1
  #  distanceNP2 = np.zeros((dims[1]**2, dims[2]-119))
   # for point in range (dims[1]**2):
        
     #   distanceNP2[point,:] = moveAverage(distanceNP[point,:], 120)
    return distanceNP
#print(distancemap(positionmap('vmd/m9-center.txt', 11))[:,30])
            
def nmfMap(distancemap,compons):   
    model = NMF(n_components = compons, max_iter=1000, tol= 1*10**-3, solver= 'mu', beta_loss= 'kullback-leibler',  )
    W = model.fit_transform(distancemap)
    H = model.components_
    print(H.shape, "H")
    print(W.shape, "W")
    for n in range (compons):
        plt.plot(W[:,n])
    plt.legend(["Component 1","Component 2"])
    plt.xlim(0,W.shape[0])
    for i in range(0, W.shape[0],int(np.sqrt(W.shape[0])) ):
                   plt.vlines(i, 0, 30, colors='Grey')
    plt.title('Protein components')
    
    plt.show()
    plt.clf()
    for indy in range(compons):
         plt.plot(H[indy,:])
         plt.title('G2 Component ' + str(indy +1))
         plt.xlabel("Frame")
         plt.show()
         plt.clf()
    np.save('tempCompons.npy', H)
    np.save('ProCompons.npy', W)
nmfMap(distancemap(positionmap('vmd/G2-center.txt', 9)  ), 2)
     

