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
   
    if reader[0] == '{':
        superlist = reader.split ('} {')
    
    #print("\n\n\n")
        superlist [0] = superlist[0].replace('{', '')
        superlist [-1] = superlist[-1].replace('}', '')
        superlist [-1] = superlist[-1][:-1]
        #print(superlist [-50:])
        #print(len(superlist))
        #print(superlist[-9:]
    else:
        superlist = reader.split('\n')
    print(len(superlist))
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
        
    # G2bonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8)]
    # G2bondsB = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7) ,(7,8), (1,6), (2,4), (1,3), (2,7)]
    # G2fbonds = [(0,1), (0,2), (2,3), (3,4), (4,5), (5,6), (3,7),  (7,8), (8,9)  ]
    # G2fbondsB = [(0,1), (0,2), (2,3), (3,4), (4,5), (5,6), (3,7),  (7,8), (8,9), (1,2), (2,4), (3,8), (3,5), (2,7)  ]
    # M9bonds = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5)]
    # M9bondsC = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5), (1,3), (3,5),  (8,10), (2,6), (1,8), (3,7), (2,9), (2,10) ]
    # M9bondsB = [ (0,1), (1,2),(2,3),(2,8), (8,9), (9,10), (3,6), (6,7),(3,4), (4,5), (1,3),(3,5), (8, 10), (2,6), (1,8), (3,7), (2,9), (2,4) ]
    # N2bonds = [(0,1), (1,2), (2,3), (3,4), (2,5), (5,6)]
    # N2bondsB = [(0,1), (1,2), (2,3), (3,4), (2,5), (5,6),(1,5), (2,4), (2,6), (1,3)]
    # N2fbonds = [(0,1), (0,2), (2,3), (3,6), (6,7), (3,4), (4,5)   ]
    # N2fbondsB = [(0,1), (0,2), (2,3), (3,6), (6,7), (3,4), (4,5), (1,2), (3,7), (2,6) , (2,4), (3,5)  ]
    # S2bonds = [(0,1), (1,2), (2,3), (3,4),(4,5), (5,6), (2,7), (7,8), (8,9), (9,10)  ]
    # S2bondsB = [(0,1), (1,2), (2,3), (3,4),(4,5), (5,6), (2,7), (7,8), (8,9), (9,10), (1,7), (8,10), (4,6), (2,4), (2,8), (1,3)  ]
    #G1M3bonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7), (6,8),(2,4), (1,6), (2,8), (1,3), (2,7) , (1,8) , (1,7), (4,6) ]
    #A looks at (2,8, 2,9, 3,9)
    #B looks at 5-11, 6-11 motion
    #G1M5bondsA =  [ (0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7), (7,8), (6,9), (9,10), (6,8), (6,10), (2,9), (2,4), (2,7), (1,3), (1,6), (4,10), (5,10), (4,9), (2,10) ]
    #G1M5bondsB =  [ (0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7), (7,8), (6,9), (9,10), (6,8), (6,10), (2,9), (2,4), (2,7), (1,3), (1,6), (2,8), (1,7), (3,8) ]
    #G2FBbonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (2,7), (7,8), (8,9), (0,10), (1,3), (2,8), (1,7), (1,10), (6,7), (1,8), (1,9), (8,10), (9,10) ]
    G3Fbonds = [(0,1), (1,2), (2,3), (3,4), (4,5), (3,6),(6,7), (2,8), (8,9), (9,10), (0,11), (1,11), (1,8), (2,4), (2,9), (1,3), (4,6), (2,6)  ]
    for bond in G3Fbonds:
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
            

    customH = np.full((numNMF,frames),.42)        


    
    return contactNP, customH, customW

def nmfMap(distancemap, customH, customW):   
    model = NMF(n_components = numNMF, max_iter=1200, tol= 1*10**-8, solver= 'mu', beta_loss= 'kullback-leibler', init ='custom' )
    W = model.fit_transform(distancemap,  W = customW, H  = customH)
    H = model.components_
    plt.plot(W)
    plt.title("Protein Binary Components")
    #plt.legend(["Component 1"])
    #plt.legend(["Component 1"]) #G1M3bonds =  [ (0,1), (1,2), (2,3), (3,4), (4,5), (2,6), (6,7), (6,8), (1,3), (2,4),(1,6), (7,8), (2,7)    ]

    plt.show()
    plt.clf()
    for indy in range(numNMF):
         plt.plot(H[indy,:], linewidth=.3)
         plt.title('G1M3 Component ' + str(indy +1))
         plt.xlabel("Frame")
         plt.show()
         plt.clf()  
    #plt.legend(["Component 1"])
    print(H.shape, "H")
    print(W.shape, "W")
    np.save('tempCompons.npy', H)
    np.save('ProCompons.npy', W)
    
        
pm =positionmap('vmd/g3f-45d.txt',12 )
nmfMap(pm[0], pm[1], pm[2])