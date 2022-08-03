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
W = np.load('ProCompons.npy')
H2 = np.zeros((1, H.shape[1]-499))
for entry in range(1):
    H2[entry,:] = moveAverage(H[entry,:], 500)
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
Hpeaks =  find_peaks(H2[0,:], prominence=.005)

print(Hpeaks)
print ("----------------\n")
print (find_peaks(H2[0,:]*-1, prominence=.005))

print ("----------------\n")
peaks = find_peaks(W[:,0], height = .7)
#print(peaks)





tupleList = []
for n in peaks[0]:
     tupleList += [ tuple( [ (n // 11) + 1, (n % 11) + 1 ])]
index = 0
for (x,y) in zip(peaks[0], peaks[1]['peak_heights']) :
        label = tupleList[index]
       #  if index %2 == 0:
       #      y = y + .5
       #      x = x -5
       # # else:
       #  #    x = x +4
        if index == 1:
             x = x-2
       #  if label == (4,2):
       #      x = x -3
        if index ==3:
            x-=1
        if index == 4 :
             x+=4
        if index == 5:
             x+=1
        if index == 6 :
             y+=0
        if index == 7 :
             x+=6

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
# peaks2 = find_peaks(W[:,1], height = 0.5)
# #print(peaks2)
# tupleList2 = []
# for n in peaks2[0]:
#       tupleList2 += [ tuple( [ (n // 11) + 1, (n % 11) + 1 ])]
# index = 0
# for (x,y) in zip(peaks2[0], peaks2[1]['peak_heights']) :
#         label = tupleList2[index]
#         if index == 3:
#             x = x +4
        
#         plt.annotate(label, # this is the text
#               (x,y), # these are the coordinates to position the label
#                   textcoords="offset points", # how to position the text
#                   xytext=(0,1), # distance from text to points (x,y)
#                   ha='center', color= "red") # horizontal alignment can be left, right or center
        
#         index += 1
        

plt.plot(W[:,0], color='green')
#plt.plot(W[:,1], color='red')
#plt.plot(W[:,1], color ='red')
#plt.ylim(0,10)

plt.title("Protein Binary Components - Component 0")
#plt.legend(['Comp 1'])

#plt.legend(peaks, tupleList)
plt.show()
plt.clf()

    
    

    
#for search in range (0,10000,10):
#search= 
#    print(H[0, search], search)
#print(search)
for n in range(1):
    plt.plot(np.linspace(0, len(H2[n,:]), len(H2[n,:])), H2[n,:])
    plt.title("S2 Component "+ str(n))
    plt.show()
    plt.clf()

for n in range(1):
    plt.plot(np.linspace(0, len(H[n,:]), len(H[n,:])), H[n,:])
    plt.title("S2 Component "+ str(n))
    plt.show()
    plt.clf()