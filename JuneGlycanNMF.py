# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:57:49 2022

@author: Abner
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

def contactReader(mapo,residu,fra):
    mapGL2 =  open(mapo, 'r')
    maptext = mapGL2.read()  
#print(maptext)
    mapGL2.close()
#print(maptext)
#print(type(maptext))
#print(type(maptext[0]))
#print(maptext.count('\n'))
    contactList = maptext.split("{")
    contactList = contactList[1:]
    print(len(contactList))
    contactNP = np.zeros((residu ** 2,fra))
    nineTimer = 0
    contactNPi = 0
    numList = []
    for stringy in contactList:
        if nineTimer  ==  residu:
            contactNP[:,contactNPi ]  = numList
            numList = []
            nineTimer = 0
            contactNPi+= 1
            
        numString = stringy.replace('}', '')
        numString = numString.replace(' ', '')
        numString = numString.replace('\n', '')
        
        #print (numString)
        for num in numString:
            numList.append(int(num))
        nineTimer += 1
    contactNP2 = contactNP.transpose()
    model = NMF(n_components=4, max_iter=500, tol= 1*10**-10, solver= 'mu',init = 'random',  beta_loss= 'kullback-leibler', alpha = .1 )#, alpha = .3  )
    W = model.fit_transform(contactNP2)
    H = model.components_
    print(H.size, "H")
    for n in range (4):
        plt.plot(H[n,:])
    plt.title('Protein components')
    plt.show()
    plt.clf()

    print ("W-size",W.shape)
    #plt.plot(W)
    #H = H.transpose()
    for indy in range(4):
         plt.plot(W[:,indy])
         plt.title('S2 Component ' + str(indy))
         plt.xlabel("Frame")
         plt.show()
         plt.clf()

    
contactReader('vmd/S2-10.txt',11,1001)

#print(len(contactNP[:,568]))    
    