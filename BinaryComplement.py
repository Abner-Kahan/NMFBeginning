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
    keepList =[(2,5), (5,2), (1,5), (5,1), (1,6), (6,1)]
    for a in range (residues):
        for b in range(residues):
            if (a,b) not in keepList:
                    zerolist.append((a,b))
                    
    
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
            

    customH = np.full((numNMF,frames),.3)        


    
    return contactNP, customH, customW

def nmfMap(distancemap, customH, customW):   
    model = NMF(n_components = numNMF, max_iter=1200, tol= 1*10**-8, solver= 'mu', beta_loss= 'kullback-leibler', init ='custom' )
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
         plt.title('N2f Component ' + str(indy +1))
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