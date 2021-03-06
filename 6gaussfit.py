#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 18:42:56 2022

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




ran1 =1600
ran2 = 1730


amide2_ran1 =1475
amide2_ran2 =1600

#1450, 1750
amide1_ran1 =1600
amide1_ran2 =1715

Humdities = [5,10,20,30,40,50,60,70,80,90,95]
solvents = ["Untreated", "MeOH", "WA45"]
broad = 15
x_range = np.linspace(ran1,ran2,ran2-ran1+1)


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
def gauss(x, h, A, x0, c):
    return h + (A * np.exp (-1*((x - x0) ** 2 / ( c ** 2))))

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

def sixgauss(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2,  \
              sigma_2, A_3, x0_3, sigma_3, A_4, x0_4, sigma_4,A_5, x0_5, sigma_5  ):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))  + \
        ( A_2 * np.exp(-(x - x0_2) ** 2 / ( sigma_2 ** 2)))   + \
    (A_3 * np.exp(-(x - x0_3) ** 2 / (sigma_3 ** 2))) +   \
       (A_4 * np.exp(-(x - x0_4) ** 2 / (sigma_4 ** 2)))  + \
           (A_5 * np.exp(-(x - x0_5) ** 2 / (sigma_5 ** 2)))
           
       
def fitFive(selec,sol):
    IRboy = fetchIr(solvents[sol]+'Sample.txt',selec,ran1,ran2)
    broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
    #IrPlotter( broadMan, f'{solvents[sol]} Spectra_old_fit', ran1,ran2)
    peak1 =1600
    peak2 =  1620 #1620  
    peak3 =  1645#1645
    peak4 =   1660 #random.randrange(ran1,ran2,1)# # random.randrange(ran1,ran2,1)
    peak5 =  1678 #random.randrange(ran1,ran2,1)
    peak6 = 1690 #1678
    print('The peaks are',peak1, peak2,peak3,peak4)
    guesses =[1] +[1,peak1,10]+[1,peak2,10]+[1,peak3,10]+[1,peak4,10]+ [1,peak5,10] + [1,peak5,10]
    constraints = ([-20,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5],[20,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15] )
                                                              
    fit_y6 = curve_fit(sixgauss,IRboy[0],IRboy[1], p0 = guesses,bounds = constraints, method ='trf')
     #print(fit_y2)
    par6 = fit_y6[0]
    yplot6 =sixgauss(x_range, par6[0],par6[1],par6[2],par6[3],par6[4],par6[5],par6[6],par6[7],par6[8],par6[9],
                     par6[10],par6[11], par6[12], par6[13],par6[14], par6[15], par6[16], par6[17], par6[18] )
     
#print("gauss2",fit_y2)
    #plt.plot(x_range,yplot6)
    #plt.xlim(max(x_range), min(x_range))
    #plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
    #plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, yplot6):.5f} ")
     
     #plt.ylim(bottom=0)
#plt.xlim(0,4000)
#fit_y2 = (curve_fit(fourgauss,IR1[0],IR1[1]))[0]
#plt.plot(IR1[0],fit_y2[1])
#plt.plot()
#plt.xlim(0,4000)
    #plt.title(f'Sum of Six Gaussians {Humdities[selec-1]}% humidity Solvent: {solvents[sol]}')

   # plt.savefig(f'Sum6Gaussians{Humdities[selec-1]}humiditySolvent{solvents[sol]}.png')
    #plt.show()
    #plt.clf()
    gauss1 = gauss(x_range, par6[0], par6[1], par6[2], par6[3])
    gauss2 = gauss(x_range, par6[0], par6[4], par6[5], par6[6])
    gauss3 = gauss(x_range, par6[0], par6[7], par6[8], par6[9])
    gauss4 = gauss(x_range, par6[0], par6[10], par6[11], par6[12])
    gauss5 = gauss(x_range, par6[0], par6[13], par6[14], par6[15])
    gauss6 = gauss(x_range, par6[0], par6[16], par6[17], par6[18])
     
    GaussOrder = np.argsort([par6[2],par6[5],par6[8],par6[11],par6[14],par6[17] ])
    gausses = [gauss1,gauss2,gauss3,gauss4,gauss5,gauss6]
    gausses =  [gausses[i] for i in GaussOrder]
     
    #plt.plot(x_range, gausses[0])
    #plt.plot(x_range, gausses[1])
    #plt.plot(x_range, gausses[2])
   # plt.plot(x_range, gausses[3])
   # plt.plot(x_range, gausses[4])
    #plt.plot(x_range, gausses[5])
    
   # plt.title(f"6 Gaussians:{Humdities[selec-1]}% humidity Solvent: {solvents[sol]}")
    #plt.legend(["Spacer", "BS","RC","AH","BT", "End Peak"])
    
    areaA = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[1], par6[2], par6[3] ) )[0]
    areaB = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[4], par6[5], par6[6] ) )[0]
    areaC = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[7], par6[8], par6[9] ) )[0]
    areaD = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[10], par6[11], par6[12] ) )[0]
    areaE = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[10], par6[11], par6[12] ) )[0]
    areaF = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[13], par6[14], par6[15] ) )[0]
    areaG = scipy.integrate.quad(gauss, ran1,ran2, (par6[0], par6[16], par6[17], par6[18] ) )[0]
    
    areas = [areaA,areaB,areaC,areaD,areaE,areaF]
    
   
    
    areas =  [areas[i] for i in GaussOrder]
    #areas = [areaA,areaB,areaC,areaD]
    
    peaks6s = ([par6[2], par6[5],par6[8],par6[11],par6[14],par6[17] ])
    peaks6s =  [peaks6s[i] for i in GaussOrder]

    
    peak6s =sorted(peaks6s)
    print (peaks6s)
    #plt.legend([f"Gauss 1 {peaks0}cm ",f"Gauss 2 {peaks1}cm ",f"Gauss 3 {peaks2}cm ",f"Gauss 4 {peaks3}cm",f"Gauss 5 {peaks4}cm"] )
    #plt.xlabel("cm^-1")
    #plt.xlim(max(x_range), min(x_range))
     #plt.ylim(bottom=0)
   # plt.savefig(f"6Gaussians{Humdities[selec-1]}humiditySolvent{solvents[sol]}.png")
   # plt.show()
    #plt.clf()
    #return peaks5s
    ind =0 
    summa =sum(areas)
    for area in areas:
        areas[ind] = area/summa
        ind +=1
    return peaks6s,areas

#for sol in range(3):
 #   for hum in range(11):
  #      fitFive(hum+1,sol)
    
#fitFive(1,0)    
#fitFive(11,0)    
#fitFive(1,1)
#fitFive(11,1)
#fitFive(1,2)
#fitFive(11,2)
indices = []
for solv in solvents:
    for humi in Humdities:
        indices.append(str(solv)+": "+str(humi)+'%')
 #       
Peaks = pd.DataFrame(columns =["Spacer", "Am I - BS","Am I -RC","Am I -AH","Am I -BT", "Am I - End Peak" ], index = indices)



model = 0
for sol in range(len(solvents)):
    for hum in range(len(Humdities)):
        
        Peaks.loc[indices[model]] = fitFive(hum+1,sol)[0]
        model+=1

Peaks.to_csv('PeaksOf6A.csv')       
        
# =============================================================================
#   
#area Finder
indices = []
for solv in solvents:
      for humi in Humdities:
          indices.append(str(solv)+": "+str(humi)+'%')
Areas = pd.DataFrame(columns =["Spacer", "Am I - BS","Am I -RC","Am I -AH","Am I -BT", "Am I - End Peak" ], index = indices)
model = 0
for sol in range(len(solvents)):
      for hum in range(len(Humdities)):          
         Areas.loc[indices[model]] = fitFive(hum+1,sol)[1]
         model+=1
# 
print(Areas)
Areas.to_csv('Areasof6A.csv')       
# =============================================================================