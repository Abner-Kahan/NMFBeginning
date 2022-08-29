#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:19:41 2022

@author: abnerkahansmack
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

numNMF =1
def moveAverage(L,n):
    return np.convolve(L, np.ones(n), 'valid') / n
print("bob")
H = np.load('tempCompons.npy')
W = np.load('ProCompons.npy')
H2 = np.zeros((numNMF, H.shape[1]-2499))
for entry in range(numNMF):
    H2[entry,:] = moveAverage(H[entry,:], 2500)
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
Hpeaks =  find_peaks(H2[0,:], prominence=.05)

#print(Hpeaks[0])
#for (p,z) in zip (Hpeaks[0], Hpeaks[1]['prominences']):
#    print (p,z)
#print ("----------------\n")
#Hmins = find_peaks(H2[0,:]*-1, prominence=1.3)
#print(Hmins[0])
#for (p,z) in zip (Hmins[0], Hmins[1]['prominences']):
 #   print (p,z)
#print ("----------------\n")
peaks = find_peaks(W[:,0], height = .2)
#print(peaks)





tupleList = []
for n in peaks[0]:
     print(n)
     tupleList += [ tuple( [ (n // 11) + 1, (n % 11) + 1 ])]
index = 0
for (x,y) in zip(peaks[0], peaks[1]['peak_heights']) :
        label = tupleList[index]

       #  if index %2 == 0:
       #      y = y + .5
       #      x = x -5
       # # else:
       #  #    x = x +4
        #if index == 1:
          #   x = x-4
       #  if label == (4,2):
       # #      x = x -3
        if index  == 0:
           x-=4
       # if index ==3:
       #    x-=5
       #  if index == 4 :
       #       x+=4
       #  if index == 5:
       #       x+=1
       #  if index == 6 :
       #       y+=0
       #  if index == 7 :
       #       x+=6

       #  if index == 8:
       #      x+=5
       #  if index == 10:
       #      x+=5
            
       # if index == 10:
       #      y  = y - 1
       # if index == 7:
        #    y = y -.5
       #     x = x - 4
        plt.annotate(label, # this is the text
             (x,y), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,3), # distance from text to points (x,y)
                 ha='center', color= "green") # horizontal alignment can be left, right or center
        index += 1

#print(peaks2)
if numNMF == 2:
    peaks2 = find_peaks(W[:,1], height = 0.5)
    tupleList2 = []
    for n in peaks2[0]:
          tupleList2 += [ tuple( [ (n // 7) + 1, (n % 7) + 1 ])]
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
        

plt.plot(W[:,0], color='green')
#plt.ylim(0,1)

#plt.plot(W[:,1], color ='red')
#plt.ylim(0,10)

plt.title("G1M3 Binary Components - Component 0")
plt.xlabel("First residue in contact with second residue")
#plt.legend(['Comp 1'])

#plt.legend(peaks, tupleList)
plt.show()
plt.clf()

    
    

    
#for search in range (0,10000,10):
#search= 
#    print(H[0, search], search)
#print(search)
for n in range(numNMF):
    plt.plot(np.linspace(0, len(H2[n,:]), len(H2[n,:])), H2[n,:], linewidth=.8)
    plt.title("G1M3 Component "+ str(n))
    plt.xlabel("Frame")
    plt.show()
    plt.clf()

for n in range(numNMF):
    plt.plot(np.linspace(0, len(H[n,:]), len(H[n,:])), H[n,:])
    plt.title("G1M3 Component "+ str(n))
    plt.show()
    plt.clf()