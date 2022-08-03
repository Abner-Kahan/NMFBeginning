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


numNMF =1

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


    
    customW = np.ones((residues**2,numNMF))      
   
    zerolist =[]
    for i in range(residues):
        zerolist.append(i*residues+i)
    G2bonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8)]
    G2bondsB = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8), (1,6), (2,4), (1,3), (2,7)]
    G2fbonds = [(0,1), (0,2), (2,3), (3,4), (4,5), (5,6), (3,7),  (7,8), (8,9)  ]
    M9bonds = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5)]
    M9bondsB = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5), (1,3),(3,5), (8, 10), (2,6), (1,8), (3,7), (2,9), (2,4) ]
    N2bonds = [(0,1), (1,2), (2,3), (3,4), (2,5), (5,6)]
    N2fbonds = [(0,1), (0,2), (2,3), (3,6), (6,7), (3,4), (4,5)   ]
    S2bonds = [(0,1), (1,2), (2,3), (3,4),(4,5), (5,6), (2,7), (7,8), (8,9), (9,10)  ]
    for bond in G2bondsB:
        zerolist.append((bond[0]*residues)+bond[1])
        zerolist.append((bond[1]*residues)+bond[0])
    print(zerolist)
    for frame in range(frames):
        for zero in zerolist:
            contactNP[zero,frame] = 0
            if frame ==0:
                customW[zero,0] =0
                if numNMF == 2 :
                    customW[zero,1] =0
                    
            
    noZero = np.count_nonzero(contactNP)
    total = contactNP.size     
    print(noZero/ total)    
            

    customH = np.full((numNMF,frames),.30)        


    
    return contactNP, customH, customW

def nmfMap(distancemap, customH, customW):   
    model = NMF(n_components = numNMF, max_iter=1200, tol= 1*10**-10, solver= 'mu', beta_loss= 'kullback-leibler', init ='custom' )
    W = model.fit_transform(distancemap,  W = customW, H  = customH)
    H = model.components_
    for n in range (numNMF):
        plt.plot(W[:,n])
    plt.title("Protein Binary Components")
    #plt.legend(["Component 1"])
    #plt.legend(["Component 1"])

    plt.show()
    plt.clf()
    for indy in range(numNMF):
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
    
        
pm =positionmap('vmd/g2-45d.txt',9 )
nmfMap(pm[0], pm[1], pm[2])