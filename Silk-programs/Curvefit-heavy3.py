#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thur Mar 10 15:16:50 2022

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


ran1 =1475
ran2 = 1705


amide2_ran1 =1475
amide2_ran2 =1600

#1450, 1750
amide1_ran1 =1600
amide1_ran2 =1710

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
       
def fourgauss(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2, sigma_2, A_3, x0_3, sigma_3):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))  + \
        ( A_2 * np.exp(-(x - x0_2) ** 2 / ( sigma_2 ** 2)))   + \
    (A_3 * np.exp(-(x - x0_3) ** 2 / (sigma_3 ** 2)))
#=============================================================================
def gaussTwo(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))
              
def gaussEight(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2, sigma_2, A_3, x0_3, sigma_3, \
                A_4, x0_4, sigma_4,  A_5, x0_5, sigma_5,  A_6, x0_6, sigma_6, A_7, x0_7, sigma_7):
    return fourgauss(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2, sigma_2, A_3, x0_3, sigma_3) +  \
            fourgauss(x, H_0, A_4, x0_4, sigma_4,  A_5, x0_5, sigma_5,  A_6, x0_6, sigma_6, A_7, x0_7, sigma_7)
def fitTwo(selec):
    IRboy = fetchIr('UntreatedSample.txt',selec,ran1,ran2)
    peak1 =  1650 #1620
    peak2 =  1550
    guesses =[0] +[1,peak1,50]+[1,peak2,50]
    constraints = ([-20,0,ran1,5,0,ran1,5],[20,np.inf,ran2,70,np.inf,ran2,70] )
    fit_y2 = curve_fit(gaussTwo,IRboy[0],IRboy[1], p0 = guesses,bounds = constraints, method ='trf')
    par2 = fit_y2[0]
    yplot2 =gaussTwo(x_range, par2[0],par2[1],par2[2],par2[3],par2[4],par2[5],par2[6] )
    plt.plot(x_range, yplot2)#,markersize=1)import matplotlib.image as mpimg
    plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
    plt.show()
    plt.clf()
    plt.plot(x_range, gauss(x_range, par2[0],par2[1],par2[2],par2[3]))
    plt.plot(x_range, gauss(x_range, par2[0],par2[4],par2[5],par2[6]))
    plt.legend(["Amide I", "Amide 2"])
#fitTwo(1)

def fitEight(selec,sol):
    IRboy = fetchIr(solvents[sol]+'Sample.txt',selec,ran1,ran2)
    broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
    IrPlotter( broadMan, f'{solvents[sol]} Spectra_old_fit', ran1,ran2)

    peak1 =  1520 #1620
    peak2 =  1545
    peak3 =  1560     # could be 1540
    peak4 =  1578
    peak5 =  1620
    peak6 =  1645
    peak7 =  1660   # could be 1640
    peak8 =  1678
    guesses =[1] +[1,peak1,10]+[1,peak2,10]+[1,peak3,10]+[1,peak4,10]+ [1,peak5,10]+ [1,peak6,10]+ [1,peak7,10]+ [1,peak8,10]
    constraints = ([-20,0,amide2_ran1, 5,0,amide2_ran1,5,0,amide2_ran1,5,0,amide2_ran1,5,0,
                    amide1_ran1,5,0,amide1_ran1,5,0,amide1_ran1,5,0,amide1_ran1,5],
                   [20,np.inf,amide2_ran2,15,np.inf,amide2_ran2,15,np.inf,amide2_ran2,15,np.inf,amide2_ran2,15,np.inf,
                    amide1_ran2,15,  np.inf, amide1_ran2,15, np.inf, amide1_ran2,15, np.inf, amide1_ran2,15] )
    fit_y8 = curve_fit(gaussEight,IRboy[0],IRboy[1], p0 = guesses,bounds = constraints, method ='trf')
    par8 = fit_y8[0]
    yplot8 =gaussEight(x_range, par8[0],par8[1],par8[2],par8[3],par8[4],par8[5],par8[6], 
                       par8[7],par8[8],par8[9],par8[10],par8[11],par8[12],par8[13],par8[14],par8[15],
                       par8[16],par8[17],par8[18],par8[19],par8[20],par8[21],par8[22],par8[23],par8[24])
    plt.plot(x_range, yplot8)#,markersize=1)import matplotlib.image as mpimg
    plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
    plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, yplot8):.5f} ")
    plt.xlim(max(x_range), min(x_range))
    plt.title(f'Sum of Eight Gaussians {Humdities[selec-1]}% humidity Solvent: {solvents[sol]}')
    plt.show()
    plt.clf()
    peaks8s = ([par8[2], par8[5],par8[8],par8[11],par8[14], par8[17], par8[20], par8[23] ])
    indexL =[]
    for peaky in sorted(peaks8s):
        indexL.append(peaks8s.index(peaky))
    amide2Colors = ['#feedde','#fdd0a2','#fdae6b','#fd8d3c','#e6550d','#a63603']
    amide1Colors = ['#f2f0f7','#dadaeb','#bcbddc','#9e9ac8','#756bb1','#54278f']
    plt.plot(x_range, (fourgauss(x_range,par8[0],par8[1],par8[2],par8[3], par8[4],par8[5],par8[6],
                                 par8[7],par8[8],par8[9], par8[10],par8[11],par8[12])), color = amide2Colors[5])
    plt.plot(x_range, gauss(x_range, par8[0],par8[1],par8[2],par8[3]),color = amide2Colors[1])
    plt.plot(x_range, gauss(x_range, par8[0],par8[4],par8[5],par8[6]),color = amide2Colors[2])
    plt.plot(x_range, gauss(x_range, par8[0],par8[7],par8[8],par8[9]),color = amide2Colors[3])
    plt.plot(x_range, gauss(x_range, par8[0],par8[10],par8[11],par8[12]),color = amide2Colors[4])

    plt.plot(x_range, (fourgauss(x_range,par8[0],par8[13],par8[14],par8[15], par8[16],par8[17],par8[18],
                                 par8[19],par8[20],par8[21], par8[22],par8[23],par8[24])),color = amide1Colors[5])
    plt.plot(x_range, gauss(x_range, par8[0],par8[13],par8[14],par8[15]),color = amide1Colors[1])
    plt.plot(x_range, gauss(x_range, par8[0],par8[16],par8[17],par8[18]),color = amide1Colors[2])
    plt.plot(x_range, gauss(x_range, par8[0],par8[19],par8[20],par8[21]),color = amide1Colors[3])
    plt.plot(x_range, gauss(x_range, par8[0],par8[22],par8[23],par8[24]),color = amide1Colors[4])


    plt.xlabel("cm^-1")
    plt.xlim(max(x_range), min(x_range))
    plt.title(f"8 Gaussians:{Humdities[selec-1]}% humidity Solvent: {solvents[sol]}")
    plt.legend(["Amide II", "BS","RC","AH","BT","Amide I","BS","RC","AH","BT"])
    
                                 
    plt.show()
    plt.clf()
    
    peaks8s = ([par8[2], par8[5],par8[8],par8[11],par8[14], par8[17], par8[20], par8[23] ])
    indexL =[]
    for peaky in sorted(peaks8s):
        indexL.append(peaks8s.index(peaky))
     
    for p in peaks8s:
        print(round(p,1), end = " ")
    


    #1620
fitEight(1,0)
#IR1 = fetchIr('WA45Sample.txt',9,ran1,ran2)
#IRBroad = gaussian_broadening(IR1,broad,ran1,ran2)
#plt.plot(IR1[0],IR1[1])
        
    