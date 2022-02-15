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
from scipy.signal import find_peaks,deconvolve

def fetchIr(path,column):
    logfile =  open(path, 'r')
    logtest = logfile.readlines()
   
    logfile.close()
    rows = len(logtest[0].split())
    columns = len (logtest)
              
    IrArray = np.zeros((rows,columns ))
    x= 0 
    y = 0
    for line in logtest:
        
        for word in line.split():
            word = word.replace(',' , '.')
            IrArray[x,y] = word
            x+=1
        y+=1
        x=0
    return IrArray[(0,column),:]

#print(fetchIr('UntreatedSample.txt',3))

def IrPlotter(item,title,ran1,ran2,leg = [], multiple = False):
    #colors = np.array([(254,230,206), (253,174,107),(230,85,13)])
    #colors = colors/256
    #colors =['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']
   # print(colors)

        
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),item,markersize=.1)
    else:
        index = 0
        for n in item:
            plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),n,markersize=.1)
            index += 1
    if len (leg) > 0:
        plt.legend(leg)
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf()    
    
ran1 =0
ran2 = 4000   
broad =20
#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
def gaussian_broadening(spectra, broaden, ran1,ran2,resolution=1):
 
    """ Performs gaussian broadening on IR spectrum
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """
    IR = np.zeros(resolution*(int(ran2-ran1) + 1))
    X = np.linspace(ran1,ran2, resolution*(int(ran2-ran1)+1))

    #for f, i in zip(spectra[:,0]):
      #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
      #  IR=np.vstack((X, IR)).T #tspec
   
    freq = spectra[0]
    inten = spectra[1]
    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
    return IR



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
Humdities = [5,10,20,30,40,50,60,70,80,90,95]
    
    
#IRF = np.zeros((33,(ran2-ran1+1)))
                   
                   
                   
#for n in range(1,12):    
  #  IRF [n-1,:] = gaussian_broadening(fetchIr('UntreatedSample.txt',n),broad,ran1,ran2)
   # IRF[n+10,:] =  gaussian_broadening(fetchIr('MeOHSample.txt',n),broad,ran1,ran2)
    #IRF[n+21,:] =  gaussian_broadening(fetchIr('WA45Sample.txt',n),broad,ran1,ran2)
IRF = np.zeros((2,(ran2-ran1+1)))
IRF [0,:] = gaussian_broadening(fetchIr('UntreatedSample.txt',3),broad,ran1,ran2)
IRF [1,:] = gaussian_broadening(fetchIr('UntreatedSample.txt',8),broad,ran1,ran2)
        
IrPlotter( IRF [0,:], 'Unstreated Spectra', ran1,ran2)  
IrPlotter( IRF [1,:], 'Unstreated Spectra', ran1,ran2)

#IrPlotter( IRF [11:22,:], 'MEOH Spectra', ran1,ran2, ['5','10','20','30','40','50','60','70','80','90','95']
   #                                                                             ,True)
#IrPlotter( IRF [22:,:], 'WA45 Spectra', ran1,ran2, ['5','10','20','30','40','50','60','70','80','90','95']
         #                                                                        ,True)