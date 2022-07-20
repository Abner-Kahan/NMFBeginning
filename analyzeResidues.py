#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:19:41 2022

@author: abnerkahansmack
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def moveAverage(L,n):
    return np.convolve(L, np.ones(n), 'valid') / n
print("bob")
H = np.load('tempCompons.npy')
W =np.load('ProCompons.npy')
H2 = np.zeros((2, H.shape[1]-499))
for entry in range(2):
    H2[entry,:] = moveAverage(H[entry,:], 500)
print(np.argmin(H2, 1))
H_sort = np.sort(H2, 1)

newList = []
for val in H_sort[0,:10]:
    place = (np.where(H2 == val))[1]
    newList += list(place)
#newarray2 = newarray
#print(newList)
#for n in range(len(newList)):
 #   for b in range(len(newList)):
        
   #     if abs(newList[n] - newList[b])<25:
   #         newList[b] = 0
#print(newList)
uniqList = []
preVal = -2000
for i in newList:
    if abs (i  - preVal) > 1000:
        uniqList.append(i)
        preVal = i
peaks = (find_peaks(W[:,1]))[0]
print(peaks)
tupleList = []
for n in peaks:
     tupleList += [ tuple( [n // 11, n % 11 ])]
print(tupleList)
    
    

    
    
    

    
    
#for search in range (0,10000,10):
#search= 
#    print(H[0, search], search)
#print(search)
for n in range(2):
    plt.plot(np.linspace(0, len(H2[0,:]), len(H2[0,:])), H2[n,:])
    plt.title("G2 Component "+ str(n))
    plt.show()
    plt.clf()