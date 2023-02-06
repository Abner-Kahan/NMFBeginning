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

def IrPlotter(item,title,ran1,ran2,ax = 0, leg = [], multiple = False):
    if ax == 0:
        plt.ylim()
    else:
        plt.ylim(top =ax)
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
    
    
#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
def gaussian_broadening(spectra, broaden, ran1,ran2,experiment =False, resolution=1):
 
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
    if experiment:
         freq = spectra[0] *.965
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



def nmf2TesterMixB(ran1,ran2,broad,res):
    #pdb.set_trace()


    Humdities = [5,10,20,30,40,50,60,70,80,90,95]
    
    
    IRF = np.zeros((33,(ran2-ran1+1)*res))
                   
                   
                   
    for n in range(1,12):
        IRF [n-1,:]= gaussian_broadening(fetchIr('UntreatedSample.txt',n),broad,ran1,ran2,res)
        IRF[n+10,:] =  gaussian_broadening(fetchIr('MeOHSample.txt',n),broad,ran1,ran2,res)
        IRF[n+21,:] =  gaussian_broadening(fetchIr('WA45Sample.txt',n),broad,ran1,ran2,res)
    
        
    
    IrPlotter( IRF [:11,:], 'Unstreated Spectra', ran1,ran2,80, ['5','10','20','30','40','50','60','70','80','90','95']
                                                                                ,True)
    IrPlotter( IRF [11:22,:], 'MEOH Spectra', ran1,ran2, 80, ['5','10','20','30','40','50','60','70','80','90','95']
                                                                                 ,True)
    IrPlotter( IRF [22:,:], 'WA45 Spectra', ran1,ran2,80, ['5','10','20','30','40','50','60','70','80','90','95']
                                                                                 ,True)
   # print(IR1.shape)
    

   # IrPlotter (gaussian_broadening(fetchIr('UntreatedSample.txt',1),broad,ran1,ran2,res), "test", ran1, ran2)
    IRF= np.transpose(IRF)
    
    model = NMF(n_components=5, max_iter=5000, tol= 1*10**-10, solver= 'mu', init='nndsvda', beta_loss= 'kullback-leibler')#, alpha = .3  )
    W = model.fit_transform(IRF)
    #IrPlotter(W[:,0])
    
    print ("W-size",W.shape)
    #print("mean", np.mean(W), np.mean(IRF[:,0]))
    numPeaks0 = (find_peaks(W[:,0])[0])+ran1
    
    numPeaks1 = (find_peaks(W[:,1])[0])+ran1
    numPeaks2 = (find_peaks(W[:,2])[0]) +ran1
    numPeaks3 = (find_peaks(W[:,3])[0]) +ran1
    numPeaks4 = (find_peaks(W[:,4])[0]) +ran1
    
    peakPlotter(W,numPeaks0,ran1,ran2,0)
    peakPlotter(W,numPeaks1,ran1,ran2,1)
    peakPlotter(W,numPeaks2,ran1,ran2,2)
    peakPlotter(W,numPeaks3,ran1,ran2,3)
    peakPlotter(W,numPeaks4,ran1,ran2,4)
    
    
    
    print ("Peaks",numPeaks0,numPeaks1,numPeaks2,numPeaks3 )
    print('nums',numPeaks0,numPeaks1,numPeaks2,numPeaks3)
    K0= np.mean([W]) /  np.mean(IRF[:,0])
    K1= np.mean([W]) /  np.mean(IRF[:,1])
    K2= np.mean([W]) /  np.mean(IRF[:,2])
    K3= np.mean([W]) /  np.mean(IRF[:,3])
    K4= np.mean([W]) /  np.mean(IRF[:,4])
   
     
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,0]*K0 ],'Output Spectra vs Untreated Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "Untreated Spectra"], True)
    
    IrPlotter(W[:,0], "Calculated 1",ran1,ran2,1.3)
    IrPlotter(W[:,1], "Calculated 2",ran1,ran2,1.3)
    IrPlotter(W[:,2], "Calculated 3",ran1,ran2,1.3)
    IrPlotter(W[:,3], "Calculated 4",ran1,ran2,1.3)
    IrPlotter(W[:,4], "Calculated 5",ran1,ran2,1.3)
    
    #IrPlotter(W[:,4], "Calculated Spectra 5",ran1,ran2)
    
    
    
    
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,1]*K1 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "MeOH A"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,2]*K2 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "MeOH B"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,3]*K3 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "WA45 A"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,4]*K4 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "WA45 B"], True)
    #IrPlotter([product[:,0],product[:,1],product[:,2],product[:,3],IRF[:,1] ],'4Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", f'MeOH Spectra 1:  {Humdities[FileSelection1]}% humidity'], True)
   # IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3] ],'4Spectra',["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra"], True)
    H = model.components_
    Hnorm = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    np.save("H.npy",H)
    np.save("Hnorm.npy",Hnorm)


nmf2TesterMixB(0,4000,20,1)  
