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

residueCount = [7,8, 9 ,10 ,11, 8, 7, 8, 8, 8, 8, 10, 11, 14, 14   ]
sugarList = ['m5-45d.txt','m6-45d.txt','m7-45d.txt', 'm8-45d.txt', 'm9-45d.txt','n1fb-45d.txt' , 
             'n2-45d.txt', 'n2b-45d.txt','n2f-45d.txt','n33-45d.txt','n36-45d.txt', 's1-45d.txt', 
             's2-45d.txt',  's33-45d.txt', 's36-45d.txt']
sugarDirectory = {'m5': 7, 'm6': 8, 'm7': 9, 'm8': 10, 'm9': 11, 'n1fb': 8, 'n2': 7, 'n2b': 8, 'n2f': 8, 'n33': 8, 'n36': 8, 's1': 10, 's2': 11, 's33': 14, 's36': 14}


numNMF = 5 #4-9
sugar = 'n36'
residues = sugarDirectory[sugar]
H = np.load('tempCompons' + '_' + sugar + '_'+ str(numNMF)+ '.npy' )
W = np.load('ProCompons' + '_' + sugar + '_'+ str(numNMF)+ '.npy' )

H_df = pd.DataFrame(data= np.transpose(H), columns =  list(range(1,numNMF+1)))
corr_matrix = H_df.corr()
Heatmap5, axH5 = plt.subplots ()
axH5 = sns.heatmap(corr_matrix,annot=True)


plt.show(axH5)
plt.clf()
def moveAverage(L,n):
    return np.convolve(L, np.ones(n), 'valid') / n

#6 is 10

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
# =================2f-10============================================================
# peaks = find_peaks(W[:,0], height = .3)
# #print(peaks)
# 
# 
# 
# 
# 
# tupleList = []Copy of N2F
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
#     ind =0 Copy of N2F
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
    


    
#for search in range (0,10000,10):Copy of N2F
#search= 
#    print(H[0, search], search)
#print(search)
for n in range(numNMF):
    plt.plot(np.linspace(0, len(H2[n,:]), len(H2[n,:])), H2[n,:], linewidth=.8)
    plt.title(  sugar + " Component "+ str(n))
    plt.xlabel("Frame")
    #plt.figure(figsize=(3.46, 2.5))
    plt.show()
    plt.clf()
    Heatmap, axH = plt.subplots ()
    axH = sns.heatmap(funTable(W[:,n],residues))
    axH.set_xlabel('Residues')
    axH.set_xticks([i+.5 for i in range(residues)] ,[i+1 for i in range(residues)])
    axH.set_ylabel('Residues')
    axH.set_yticks([i + .5 for i in range(residues)] ,[i+1 for i in range(residues)])
    #plt.figure(figsize=(3.46, 2.5))
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
    
def getPercentages(L):
    Sentence = ''
    for we in range(len(L)):
        if L[we] != 0 :
            Sentence += f" {L[we]} % Component {we}   "
    return Sentence
for n in range(numNMF):
    print(np.max(H[n]/ np.sum(H,axis=0)))
print("__________________________________")
for n in range(numNMF):
    q = np.argmax(H[n]/ np.sum(H,axis=0))
    print (q)
    print(getPercentages(np.round(( H[:,q] / np.sum(H[:,q])*100),2)))
    print('\n')