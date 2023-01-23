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

mapGL2 =  open('10incr.txt', 'r')
maptext = mapGL2.read()  
#print(maptext)
mapGL2.close()
#print(maptext)
#print(type(maptext))
#print(type(maptext[0]))
#print(maptext.count('\n'))
contactList = maptext.split("{")
contactList = contactList[1:]
#print(len(contactList))
contactNP = np.zeros((81,1000))
nineTimer = 0
contactNPi = 0
numList = []
for stringy in contactList:
    if nineTimer  ==  9:
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
model = NMF(n_components=3, max_iter=500, tol= 1*10**-10, solver= 'mu',init = 'random',  beta_loss= 'kullback-leibler', alpha = .1 )#, alpha = .3  )
W = model.fit_transform(contactNP2)
H = model.components_

print ("W-size",W.shape)
#plt.plot(W)
#H = H.transpose()
for indy in range(3):
     plt.plot(W[:,indy])
     plt.title('G2 Component ' + str(indy))
     plt.xlabel("Frame")
     plt.show()
     plt.clf()

    


#print(len(contactNP[:,568]))    
    