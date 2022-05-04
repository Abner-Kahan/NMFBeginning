# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 08:20:12 2022

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


amide1_ran1 =1600
amide1_ran2 =1715


Humdities = [5,10,20,30,40,50,60,70,80,90,95]
solvents = ["Untreated", "MeOH", "WA45"]
broad = 15
x_range = np.linspace(amide1_ran1,amide1_ran2,(amide1_ran2-amide1_ran1+1))

def fetchIr(path):
    logfile =  open(path, 'r')
    logtest = logfile.readlines()
   
    logfile.close()
    
    rows = len(logtest)
    columns = 2
              
    IrArray = np.zeros((rows,columns ))
    
    x = 0
    for line in logtest:
        line2 = line.split()
        IrArray[x] = line2
        x +=1

    
    
    
    return IrArray


#fetch
#print(fetchIr('UntreatedSample.txt',3))
print(fetchIr('BGH2c_h.dat'))
def IrPlotter(item,title,ran1,ran2, leg = [], multiple = False):
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,len(item)),item,markersize=.1)
    else:
        for n in item:
            plt.plot(np.linspace(ran1,ran2,len(n)),n,markersize=.1)
    if len (leg) > 0:
        plt.legend(leg,fontsize='x-small')
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")    
    plt.show()
    plt.clf() 

#plt.zeros(4000)
scal = np.average(fetchIr('BGH2c_h.dat')[:,1]) / np.average(fetchIr('BGH2h_h.dat')[:,1])
plt.plot(fetchIr('BGH2c_h.dat')[:,0],fetchIr('BGH2c_h.dat')[:,1])
plt.plot(fetchIr('BGH2h_h.dat')[:,0],fetchIr('BGH2h_h.dat')[:,1]*scal)
plt.xlabel("cm^-1")    
plt.ylim(0,1.1)
plt.xlim(1575,1750)
plt.title('Scaled BGH2 HighFreq')


plt.legend(["cold", "hot"])
plt.show()
plt.clf() 

      
         
         
plt.plot(fetchIr('BGH2c_l.dat')[:,0],fetchIr('BGH2c_l.dat')[:,1])
plt.title('BGH2 LowFreq')

plt.xlabel("cm^-1")    
plt.ylim(0,1.1)
plt.xlim(1000,1200)

plt.plot(fetchIr('BGH2h_l.dat')[:,0],fetchIr('BGH2h_l.dat')[:,1])

plt.legend(["cold", "hot"])
plt.show()
plt.clf() 

indList = []

#subtable  = np.linspace(min(fetchIr('BGH2c_h.dat')[:,0]), max(fetchIr('BGH2c_h.dat')[:,0]), 2000)
# for i in fetchIr('BGH2c_h.dat')[:,0] :
#     print(i)
#     mini = 10
#     ind = 0
#     inde = -1
#     for b in fetchIr('BGH2h_h.dat')[:,0]:
#         print(b)
#         if abs(i-b) < mini:
#              mini = abs(i-b)
#              inde = ind
#              print(inde)
#         ind +=1
#indList += [inde]
#print(inde)
#x1 = fetchIr('BGH2c_h.dat')[:,0]
#x2 = fetchIr('BGH2h_h.dat')[:,0]
#y2=  fetchIr('BGH2h_h.dat')[:,1]
#y1= fetchIr('BGH2c_h.dat')[:,1]
#fetchIr('BGH2h_h.dat')[:,1]

#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
# def gaussian_broadening(spectra, broaden, experiment =False, resolution=1):
 
#     """ Performs gaussian broadening on IR spectrum
#     generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
#     spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
#     :param broaden: (float) gaussian broadening in wn-1
#     :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
#     """
#     ran2= np.ceil(np.max(spectra))
#     ran1= np.floor(np.min(spectra[:,0]))
#     IR = np.zeros(resolution*(int(ran2-ran1) + 1))
#     X = np.linspace(ran1,ran2, resolution*(int(ran2-ran1)+1))

#     #for f, i in zip(spectra[:,0]):
#       #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
#       #  IR=np.vstack((X, IR)).T #tspec
   
#     freq = spectra[:,0]
#     if experiment:
#          freq = spectra[0] *.965
#     inten = spectra[:,1]
#     #print(len(freq))
#     for f,i in zip(freq,inten):
#        IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
#     return IR,ran1,ran2
# specFull = gaussian_broadening(fetchIr('BGH2c_h.dat'), 1)

#IrPlotter(specFull[0], "Plot 1", specFull[1], specFull[2])
#Untreated =gaussian_broadening(fetchIr('UntreatedSample.txt',8),25,1000,4000)
#IrPlotter(Untreated,"untreated",1000,4000)

#print(isinstance(fileList[0],str))
#plt.vlines(file2Spectra(fileList[1])[:,0],0,file2Spectra(fileList[0])[:,1])

def peakPlotter(W,numPeaks,ran1,ran2,col):
    plt.plot(np.linspace(ran1,ran2,(ran2-ran1+1)), W[:,col])
    plt.vlines(numPeaks, 0,W[:,col][numPeaks-ran1])
    plt.legend(np.ndarray.tolist(numPeaks))
    plt.ylim(top =2)
    plt.title(f"peaks Calculated Spectra{col}")
    plt.show()
    plt.clf()
def SortNp(npArray):
    Difss = [tuple(x) for x in npArray]
    Difss = sorted(Difss)
#DifPlot = sorted(DifPlot)
    DiffList = []
    for t in Difss:
        DiffList.append(t)
    DiffList = np.array(DiffList)
    return DiffList
def nmf2TesterMixB(coldSpec, HotSpec):
    coldSpec = fetchIr(coldSpec)
    
    HotSpec = fetchIr(HotSpec)
    coldX= coldSpec[:,0]
    hotX = HotSpec[:,0]
    coldY = coldSpec[:,1]
    hotY = HotSpec[:,1]
    #HotYguess = np.interp(np.linspace(amide1_ran1,amide1_ran2,116),coldX,coldY)
    HotYguess = np.interp(hotX,coldX,coldY)
    #HotYguess is the predicted values of the cold spectra at the Hot's X Values
    #ColdYguess = np.interp(np.linspace(amide1_ran1,amide1_ran2,116),hotX,hotY)
    ColdYguess = np.interp(coldX,hotX,hotY)
     #ColdYguess is the predicted values of the hot spectra at the Cold's X Values
    print("Hot prediction lenghth", len(HotYguess), "Cold Pediction length", len(ColdYguess),
          "Cold", len(coldY), "Hot", len(hotY))
    HotGuess = np.stack((hotX,HotYguess),axis=-1)
    #print(HotGuess)
   
    
    ColdGuess = np.stack((coldX,ColdYguess),axis=-1)
    HotWhole=np.concatenate((HotSpec, ColdGuess), axis = 0)
    ColdWhole=np.concatenate((coldSpec, HotGuess), axis = 0)
    #print(HotWhole.shape, ColdWhole.shape)
    ColdWhole = SortNp(ColdWhole)
    HotWhole = SortNp(HotWhole)
    #
    #plt.plot(ColdWhole[:,0],ColdWhole[:,1])
    #plt.plot(HotWhole[:,0],HotWhole[:,1])
    
    
    #plt.title("Input plots")
    #plt.legend(["Cold", "Hot"])
    #plt.xlabel("cm^-1")    
#
    #plt.show()
    #plt.clf()
    #IRF = np.stack((HotYguess,ColdYguess),axis=-1)
    IRF = np.stack((HotWhole[:,1], ColdWhole[:,1]), axis = -1)
    print(IRF.shape)
    #plt.plot(specFull[:,0], specFull[:,1])
    
    # = np.array(spec1,spec2)
    
    #IRF= np.transpose(IRF)
    # print(IRF.shape)
    model = NMF(n_components=2, max_iter=1000, tol= 1*10**-10, solver= 'mu',  beta_loss= 'kullback-leibler', init= 'random')#, alpha = .3  )
    W = model.fit_transform(IRF)
    np.savetxt("NMFHi.csv", W, delimiter=",")
    H = model.components_
    Reoncs = (np.matmul(W,H).shape)
    
    
    plt.plot(W[:,0])
    plt.xlabel("cm^-1")    
    plt.title("NMF highFreq Output plots 1/2 ")
    plt.show()
    plt.clf()
    plt.plot(W[:,1])
    plt.xlabel("cm^-1")    
    plt.title("NMF highFreq Output plots 2/2 ")
    plt.show()
    plt.clf()
    # plt.plot(np.linspace(amide1_ran1,amide1_ran2,116),W[:,2])
    # plt.xlabel("cm^-1")    
    # plt.title("NMF highFreq Output plots 3/4 Low")
    # plt.show()
    # plt.clf()
    # plt.plot(np.linspace(amide1_ran1,amide1_ran2,116),W[:,3])

    # plt.xlabel("cm^-1")    
    # plt.title("NMF highFreq Output plots 4/4 Low")
    
    # plt.show()
    # plt.clf()
    print(H.shape)
    print(W.shape)
    # #plt.vlines(W[0], W[1])
    # print(H , '\n\n')
    # print(W, '\n\n')
    # print(len(np.matmul(W,H)[0]))
    # plt.plot(np.matmul(W,H)[0],  np.matmul(W,H)[1] )
nmf2TesterMixB('BGH2c_h.dat', 'BGH2h_h.dat')
#sumsol  = np.load('sumsol.npy')
