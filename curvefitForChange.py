#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:16:50 2022

@author: abnerkahan
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scipy
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import savgol_filter
from scipy.signal import find_peaks,deconvolve, peak_widths
from scipy.optimize import curve_fit
import matplotlib.image as mpimg
import warnings


ran1 =1600
ran2 = 1700
broad = 15


def ImaginaryEnergy(spectra):
    peaks, _ = find_peaks(spectra)
    results_half = peak_widths(spectra, peaks, rel_height=0.5)
    ImaginaryEnergy = np.average (results_half[0])/2
    return ImaginaryEnergy
def Y(wavelength, spectra,V_oi):
    L = scipy.misc.derivative (lambda x:(spectra[int(x)]), wavelength)/spectra[int(wavelength)]
    return (L**-1)/(L**-2 + V_oi**2)

def num(wavelength, TheoSpec,ExperSpec,V_oi):
    return (Y(wavelength, TheoSpec,V_oi) - Y(wavelength, ExperSpec,V_oi))**2
def denom(wavelength, TheoSpec,ExperSpec,V_oi):
    return (Y(wavelength, TheoSpec,V_oi)**2) + (Y(wavelength, ExperSpec,V_oi)**2)

def ypendry(TheoSpec,ExperSpec):
    #specDif = SpecDifferenceSq(TheoSpec,ExperSpec)
    return ( ( scipy.integrate.quad(num,0, len(TheoSpec)-1, (TheoSpec,ExperSpec, ImaginaryEnergy(TheoSpec)), limit=100)[0]) / (
    scipy.integrate.quad(denom,0, len(TheoSpec)-1, (TheoSpec,ExperSpec,ImaginaryEnergy(TheoSpec)), limit=100))[0])



def fetchIr(path,column,ran1,ran2):
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
            if (x == 0) and (float(word) < ran1 or float(word) > ran2):
                
                break


            #print(word, '\n')
            #if float(word) < ran1:
          #      break
           # if float(word) > ran2:
            #    break
            IrArray[x,y] = word
            x+=1
        y+=1
        x=0
    #print(IrArray.shape)
    mask = (IrArray[0,:]) != [0] *len(IrArray[0])
    #print(mask)
    
    IrArray2 = IrArray[:,mask]
    #print(IrArray2.shape)
    #plt.scatter(IrArray2[0],IrArray2[1])
    #print(IrArray)
    
    
    return IrArray2[(0,column),:]

def gaussian_broadening(spectra, broaden, ran1,ran2,resolution=1,theory=False):
 
    """ Performs gaussian broadening on IR spectrumimport matplotlib.image as mpimg
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


def IrPlotter(item,title,ran1,ran2,leg = [], multiple = False):
    #colors = np.array([(254,230,206), (253,174,107),(230,85,13)])
    #colors = colors/256
    #colors =['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']
   # print(colors)

        
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),item)
    else:
        index = 0
        for n in item:
            plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),n)
            index += 1
    if len (leg) > 0:
        plt.legend(leg)
    #plt.ylim(bottom=0)
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.xlim(ran2,ran1)
    
    plt.show()
    plt.clf()   

def gauss(x, h, A, x0, c):
    return h + (A * np.exp (-1*((x - x0) ** 2 / ( c ** 2))))

#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
    
def fourgauss(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2, sigma_2, A_3, x0_3, sigma_3):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))  + \
        ( A_2 * np.exp(-(x - x0_2) ** 2 / ( sigma_2 ** 2)))   + \
       (A_3 * np.exp(-(x - x0_3) ** 2 / (sigma_3 ** 2)))
 
b=11
Humdities = [5,10,20,30,40,50,60,70,80,90,95]
for n in range (1):
      IRboy = fetchIr('UntreatedSample.txt',b,ran1,ran2)
      broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
      print("Run", n, "\n\n\n\n\n")
      peak1 =  random.randrange(ran1,ran2,1) #1620
      peak2 =  random.randrange(ran1,ran2,1) #1645
      peak3 =   1660# # random.randrange(ran1,ran2,1)
      peak4 =  random.randrange(ran1,ran2,1) #1678
      print('The peaks are',peak1, peak2,peak3,peak4)
      guesses =[1] +[1,peak1,10]+[1,peak2,10]+[1,peak3,10]+[1,peak4,10]
      constraints = ([-20,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5],[20,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15] )
      x_range = np.linspace(ran1,ran2,ran2-ran1+1)                                                         
      fit_y2 = curve_fit(fourgauss,IRboy[0],IRboy[1], p0 = guesses,bounds = constraints, method ='trf')
      #print(fit_y2)
      par4 = fit_y2[0]
      #erro = par4[1]
      #print(fit_y2[1])
      minSpec = min(fourgauss(x_range, par4[0],par4[1],par4[2],par4[3],par4[4],par4[5],par4[6],par4[7],par4[8],par4[9],par4[10],par4[11], par4[12] ))
      max1 = max (gauss(x_range, par4[0], par4[1], par4[2], par4[3]))
      max2 = max(gauss(x_range, par4[0], par4[4], par4[5], par4[6]))
      max3 = max (gauss(x_range, par4[0], par4[7], par4[8], par4[9]))
      max4 = max(gauss(x_range, par4[0], par4[10], par4[11], par4[12]))
      br = False
      for maxNum in [max1,max2,max3,max4]:
          #print("XD \n", maxNum / (min([max1,max2,max3,max4])-minSpec))
          if maxNum / abs((min([max1,max2,max3,max4])-minSpec)) > 10:
              br = True
      if br:
          print('bad fit')#continue
      if minSpec < 0:
          continue
         
     
#print("fit2",fit_y2)
#print("bob", par4)
      yplot4 =fourgauss(x_range, par4[0],par4[1],par4[2],par4[3],par4[4],par4[5],par4[6],par4[7],par4[8],par4[9],par4[10],par4[11], par4[12] )
     
#print("gauss2",fit_y2)
      plt.plot(x_range,yplot4)
      plt.xlim(max(x_range), min(x_range))
      plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
      plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, yplot4):.5f} ")
     
      #plt.ylim(bottom=0)
#plt.xlim(0,4000)
#fit_y2 = (curve_fit(fourgauss,IR1[0],IR1[1]))[0]
#plt.plot(IR1[0],fit_y2[1])
#plt.plot()
#plt.xlim(0,4000)
      plt.title(f'RUN {n}: Sum of Four Gaussians {Humdities[b-1]}% humidity')
     
      plt.show()
      plt.clf()
      #print (max (gauss(x_range, par4[0], par4[1], par4[2], par4[3])))
    # print (max (gauss(x_range, par4[0], par4[4], par4[5], par4[6])))
    # print (max (gauss(x_range, par4[0], par4[7], par4[8], par4[9])))
    # print (max (gauss(x_range, par4[0], par4[10], par4[11], par4[12])))

#OneOf4_gauss = 

     
      peaks4s = ([par4[2], par4[5],par4[8],par4[11]])
      indexL =[]
      peaks4s = sorted(peaks4s)
    
     
      print(peaks4s)
     
                                 
      gaussA = gauss(x_range, par4[0], par4[1], par4[2], par4[3])
      #area1 = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[1], par4[2], par4[3] ) )[0]
      
      gaussB = gauss(x_range, par4[0], par4[4], par4[5], par4[6])
      gaussC = gauss(x_range, par4[0], par4[7], par4[8], par4[9])
      gaussD = gauss(x_range, par4[0], par4[10], par4[11], par4[12])
      
      areaA = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[1], par4[2], par4[3] ) )[0]
      areaB = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[4], par4[5], par4[6] ) )[0]
      areaC = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[7], par4[8], par4[9] ) )[0]
      areaD = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[10], par4[11], par4[12] ) )[0]
      areas = [areaA,areaB,areaC,areaD]
      
      gausses = [gaussA,gaussB,gaussC,gaussD]
      GaussOrder = np.argsort([par4[2],par4[5],par4[8],par4[11]])
      
      gausses =  [gausses[i] for i in GaussOrder]
      areas =  [areas[i] for i in GaussOrder]
      
      #while len (gaussCopy>1):
        #  maxList += gaussCopy.index(max(gaussCopy))
          
    
                            

             
      plt.plot(x_range, gausses[0])
      plt.plot(x_range, gausses[1])
      plt.plot(x_range, gausses[2])
      plt.plot(x_range, gausses[3])                            
         
     
         
     
      plt.title(f"Four Gaussians: ({peaks4s[0]:4.1f}, {peaks4s[1]:4.1f}, {peaks4s[2]:4.1f}, {peaks4s[3]:4.1f} {Humdities[b-1]}% humidity)")
      #plt.legend(["Beta Sheet", "Random Coil","Alpha Helix","Beta Turn"])
      plt.legend([f"Gauss1 {areas[0]:.4}", f"Gauss2 {areas[1]:.4}", f"Gauss3 {areas[2]:.4}", f"Gauss4 {areas[3]:.4}"])
      plt.xlabel("cm^-1")
      plt.xlim(max(x_range), min(x_range))
      #plt.ylim(bottom=0)
      plt.show()
      plt.clf()
      #IRF [0,:] = gaussian_broadening(IR1,broad,ran1,ran2)
      #IRF [1,:] = gaussian_broadening(IR2,broad,ran1,ran2)   
      print(f"Error: {ypendry(broadMan, yplot4):.5f}")
      print ("Areas", areas, "\n\n\n\n")
