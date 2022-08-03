#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 09:45:28 2022

@author: abnerkahansmack
"""
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import NMF
from scipy.spatial.distance import cdist

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


            
            
    
    zerolist =[]
    for i in range(residues):
        zerolist.append(i*residues+i)
    #G2bonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8)]
    #G2fbonds = [(0,1), (0,2), (2,3), (3,4), (4,5), (5,6), (3,7),  (7,8), (8,9)  ]
    M9bonds = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5)]
    M9bondsB = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5), (1,3),(3,5), (8, 10) ]
    #N2bonds = [(0,1), (1,2), (2,3), (3,4), (2,5), (5,6)]\
   # N2fbonds = [(0,1), (0,2), (2,3), (3,6), (6,7), (3,4), (4,5)   ]
    #S2bonds = [(0,1), (1,2), (2,3), (3,4),(4,5), (5,6), (2,7), (7,8), (8,9), (9,10)  ]
    for bond in M9bondsB:
        zerolist.append((bond[0]*residues)+bond[1])
        zerolist.append((bond[1]*residues)+bond[0])
    print(zerolist)
    for frame in range(frames):
        for zero in zerolist:
            contactNP[zero,frame] = 0
    
    return contactNP
def nmfMap(distancemap):   
    model = NMF(n_components =2, max_iter=1200, tol= 1*10**-10, solver= 'mu', beta_loss= 'kullback-leibler', init ='random')
    W = model.fit_transform(distancemap)
    H = model.components_
    for n in range (2):
        plt.plot(W[:,n])
    plt.title("Protein Binary Components")
    #plt.legend(["Component 1"])
    #plt.legend(["Component 1"])

    plt.show()
    plt.clf()
    for indy in range(2):
         plt.plot(H[indy,:], linewidth=.3)
         plt.title('S2 Component ' + str(indy +1))
         plt.xlabel("Frame")
         plt.show()
         plt.clf()  
    #plt.legend(["Component 1"])
    print(H.shape, "H")
    print(W.shape, "W")
    np.save('tempCompons.npy', H)
    np.save('ProCompons.npy', W)
    
        
nmfMap(positionmap('vmd/m9-45d.txt',11 ))