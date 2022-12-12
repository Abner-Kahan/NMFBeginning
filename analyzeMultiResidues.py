#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:19:41 2022

@author: abnerkahansmack
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import pandas as pd
import seaborn as sns


numNMF =6
residues = 11
def moveAverage(L,n):
    return np.convolve(L, np.ones(n), 'valid') / n
print("bob")
H = np.load('tempCompons.npy')
W = np.load('ProCompons.npy')
H2 = np.zeros((numNMF, H.shape[1]-999))
for entry in range(numNMF):
    H2[entry,:] = moveAverage(H[entry,:], 1000)
#print(np.argmin(H2, 1))
#H_sort = np.sort(H2, 1)

#newList = []
#for val in H_sort[0,:10]:
#    place = (np.where(H2 == val))[1]
#    newList += list(place)
#newarray2 = newarray
#print(newList)
#for n in range(len(newList)):
 #   for b in range(len(newList)):
        
   #     if abs(newList[n] - newList[b])<25:
   #         newList[b] = 0
#print(newList)
#uniqList = []
#preVal = -2000
#for i in H_sort[0, :]:
 #   if abs (i  - preVal) > 1000:
 #       uniqList.append(i)
 #       preVal = i
Hpeaks =  find_peaks(H2[0,:], prominence=.01)
#print(Hpeaks[0])
#for (p,z) in zip (Hpeaks[0], Hpeaks[1]['prominences']):
#    print (p,z)
#print ("----------------\n")
#Hmins = find_peaks(H2[0,:]*-1, prominence=1.3)
#print(Hmins[0])
#for (p,z) M6Bondsin zip (Hmins[0], Hmins[1]['prominences']):
 #   print (p,z)
#print ("----------------\n")
# =============================================================================
# peaks = find_peaks(W[:,0], height = .3)
# #print(peaks)
# 
# 
# 
# 
# 
# tupleList = []
# for n in peaks[0]:
#      print(n)
#      tupleList += [ tuple( [ (n // 11 + 1), (n % 11 + 1 )])]
# index = 0
# for (x,y) in zip(peaks[0], peaks[1]['peak_heights']) :
#         label = tupleList[index]
# 
# 
#         plt.annotate(label, # this is the text
#              (x,y), # these are the coordinates to position the label
#                  textcoords="offset points", # how to position the text
#                  xytext=(0,3), # distance from text to points (x,y)
#                  ha='center', color= "green") # horizontal alignment can be left, right or center
#         index += 1
# =============================================================================

#print(peaks2)
if numNMF == 2:
    peaks2 = find_peaks(W[:,1], height = 0.5)
    tupleList2 = []
    for n in peaks2[0]:
          tupleList2 += [ tuple( [ (n // 11) + 1, (n % 11) + 1 ])]
    index = 0
    for (x,y) in zip(peaks2[0], peaks2[1]['peak_heights']):
            label = tupleList2[index]
            if index == 3:
                x = x +4
            
            plt.annotate(label, # this is the text
                  (x,y), # these are the coordinates to position the label
                      textcoords="offset points", # how to position the text
                      xytext=(0,1), # distance from text to points (x,y)
                      ha='center', color= "red") # horizontal alignment can be left, right or center
            
            index += 1
            
    plt.plot(W[:,1], color='red')
        

#plt.plot(np.linspace(0,11**2,121), W, color='green')
#plt.ylim(0,1)

#plt.plot(W[:,1], color ='red')
#plt.ylim(0,10)

#plt.title("M9 Binary Components - Component 0")
#plt.xlabel("First residue in contact with second residue")
#plt.legend(['Comp 1'])

#plt.legend(peaks, tupleList)
#plt.show()
#plt.clf()

def funTable(W2,resids):
# =============================================================================
#     ext = np.zeros((resids,resids))
#     ind =0 
#     for a in range(resids):
#         for b in range(resids):
#             ext[a,b]= W2[ind]
#             ind+=1
#     print(ext)
#     return plt.table(cellText = ext)
# =============================================================================
    mapContact = pd.DataFrame(data = [W2[:resids]])
    #print(mapContact)
# =============================================================================
    for n in range(resids-1):
         n=n+1
         mapContact.loc[n] = W2[n*resids:(n+1)*resids]
    return mapContact
    
# =============================================================================
    


    
#for search in range (0,10000,10):
#search= 
#    print(H[0, search], search)
#print(search)
for n in range(numNMF):
    plt.plot(np.linspace(0, len(H2[n,:]), len(H2[n,:])), H2[n,:], linewidth=.8)
    plt.title("G1M5 Component "+ str(n))
    plt.xlabel("Frame")
    plt.show()
    plt.clf()
    Heatmap, axH = plt.subplots ()
    axH = sns.heatmap(funTable(W[:,n],residues))
    axH.set_xlabel('Residues')
    axH.set_ylabel('Residues')
    plt.show()
    plt.clf()

# =============================================================================
#     peaks = find_peaks(W[:,n], height = 0.5)
#     tupleList =[]
#     for i in peaks[0]:
#         tupleList += [ tuple( [ (i // 11 + 1), (i % 11 + 1 )])]
#     index = 0
#     for (x,y) in zip(peaks[0], peaks[1]['peak_heights']) :
#             label = tupleList[index]
#             plt.annotate(label, # this is the text
#                  (x,y), # these are the coordinates to position the label
#                      textcoords="offset points", # how to position the text
#                      xytext=(0,3), # distance from text to points (x,y)
#                      ha='center', color= "green") # horizontal alignment can be left, right or center
#             index += 1
# =============================================================================

    plt.show()
    plt.clf()
    
    
#for n in range (50):
 #   print(.02+n/100, np.where(abs(H - .02-n/100) <.01))
 
#numpy.random.choice
def myRandom(npAr):
    if npAr.size > 0:
        return np.random.choice(npAr)
    else:
        return ""
    
#for n in range (50):
 #   print(0+n/100, myRandom(np.where(abs(H - 0-n/100) <.01)[1]))


# for n in [7947,92768,85962, 64950 ,53011 ,83369 ,2078 , 45832 ,55691 ,57406 ,56477]:
#     print(f"H[{n}] = {H[0, n]}")
 
 #np.where(abs(H[0,20000:]-.38)<.01)
# =============================================================================
# counts = []
# for n in range (int(np.ceil(np.max(H))*10)+1):
#     counts.append(np.count_nonzero(np.logical_and (H[0] <= n/10, H[0] > (n-1)/10 )))
# print(counts)
# plt.bar(np.linspace(0,3, 31), counts, log = True)
# plt.xlabel('value')
# plt.ylabel('Frequency')
#      
# =============================================================================
