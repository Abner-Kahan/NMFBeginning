#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:16:50 2022

@author: abnerkahan
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
from scipy.optimize import curve_fit


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

ran1 =1500
ran2 = 1750   
broad =20
#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
def gaussian_broadening(spectra, broaden, ran1,ran2,resolution=1,theory=False):
 
    """ Performs gaussian broadening on IR spectrum
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1"""
    """
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """
    IR = np.zeros(resolution*(int(ran2-ran1) + 1))
    X = np.linspace(ran1,ran2, resolution*(int(ran2-ran1)+1))

    #for f, i in zip(spectra[:,0]):
      #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
      #  IR=np.vstack((X, IR)).T #tspec
   
    freq = spectra[0]
    inten = spectra[1]
    if theory:
        freq = spectra[:,0]*.965
        inten = spectra[:,1]



    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
    return IR
#Define the Gaussian function


IRF = np.zeros((2,(ran2-ran1+1)))
IR1 = fetchIr('UntreatedSample.txt',3)
IR2 = fetchIr('UntreatedSample.txt',3)

def gauss(x, A, x0):
    return A * np.exp(-(x - x0) ** 2 / (2 * 20 ** 2))
def fourgauss(x, H_0, A_0, x0_0, sigma_0, x_0, H_1, A_1, x0_1, sigma_1, x_1, H_2, A_2, x0_2, sigma_2, x_2, H_3, A_3, x0_3, sigma_3):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (2 * sigma_0 ** 2))) +   \
             (H_1 + A_1 * np.exp(-(x - x0_1) ** 2 / (2 * sigma_1 ** 2)))  + \
        (H_2 + A_2 * np.exp(-(x - x0_2) ** 2 / (2 * sigma_2 ** 2)))   + \
    ( H_3 + A_3 * np.exp(-(x - x0_3) ** 2 / (2 * sigma_3 ** 2)))

fit_y = curve_fit(gauss,IR1[0],IR1[1],method='dogbox')
print("fit1", fit_y[0])
#fit_y =gauss(IR1[0], fit_y[0],fit_y[1],fit_y[2],fit_y[3])\
y_out =[]
#for x in IR1[0]:
 #   y_out.append(gauss(x,fit_y[0],fit_y[1],fit_y[2],fit_y[3]))
#plt.plot(IR1[0],y_out)
plt.plot(IR1[0], gauss(IR1, *fit_y[0])[0])
#plt.xlim(0,4000)
plt.title("fit")
plt.show()
plt.clf()


fit_y2 = (curve_fit(fourgauss,IR1[0],IR1[1]))[0]
print("fit2", fit_y2)
fit_y2 =fourgauss(IR1, fit_y2[0],fit_y2[1],fit_y2[2],fit_y2[3],fit_y2[4],fit_y2[5],fit_y2[6],fit_y2[7],fit_y2[8],fit_y2[9],
                  fit_y2[10],fit_y2[11],fit_y2[12],fit_y2[13],fit_y2[14],fit_y2[15],fit_y2[16],fit_y2[17], fit_y2[18] )
print("gauss2",fit_y2)
plt.plot(IR1[0],fit_y2[1])
plt.plot()
#plt.xlim(0,4000)
plt.title("fit2")
plt.show()
plt.clf()
IRF [0,:] = gaussian_broadening(IR1,broad,ran1,ran2)
IRF [1,:] = gaussian_broadening(IR2,broad,ran1,ran2)   

IrPlotter( IRF [0,:], 'Unstreated Spectra_Old_fit', ran1,ran2)  
IrPlotter( IRF [1,:], 'Unstreated Spectra_new_fit', ran1,ran2)