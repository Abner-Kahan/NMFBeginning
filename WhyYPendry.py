#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 19:06:48 2022

@author: abnerkahansmack
"""
import time
import glob
import warnings
#warnings.filterwarnings("ignore", message="/home/abnerkahan/anaconda3/lib/python3.9/site-packages/sklearn/decomposition/_nmf.py:1637: ConvergenceWarning:e")
#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF
from scipy.signal import chirp, find_peaks, peak_widths
import scipy
import re

timeA = time.perf_counter_ns()
def ImaginaryEnergy(spectra):
    peaks, _ = find_peaks(spectra)
    results_half = peak_widths(spectra, peaks, rel_height=0.5)
    ImaginaryEnergy = np.average (results_half[0])/2
    #print(ImaginaryEnergy)
    return ImaginaryEnergy

def Y(wavelength, spectra,V_oi):
    L = scipy.misc.derivative (lambda x:(spectra[int(x)]), wavelength)/spectra[int(wavelength)]
    #if L == 0:
        #print(str(wavelength))
    return (L**-1)/(L**-2 + V_oi**2)
def newY(wavelength, spectra,V_oi):
    L = (spectra[wavelength+1]-spectra[wavelength])/spectra[wavelength]
    return (L**-1)/(L**-2 + V_oi**2)
def num(wavelength, TheoSpec,ExperSpec,V_oi):
    return (Y(wavelength, TheoSpec,V_oi) - Y(wavelength, ExperSpec,V_oi))**2
def denom(wavelength, TheoSpec,ExperSpec,V_oi):
    return (Y(wavelength, TheoSpec,V_oi)**2) + (Y(wavelength, ExperSpec,V_oi)**2)

def ypendry(TheoSpec,ExperSpec):
    eTheo = ImaginaryEnergy(TheoSpec)
    #eExper = ImaginaryEnergy(ExperSpec)
    #specDif = SpecDifferenceSq(TheoSpec,ExperSpec)
    return ( ( scipy.integrate.quad(num,0, len(TheoSpec)-1, (TheoSpec,ExperSpec, eTheo),limit =100)[0]) / (
    scipy.integrate.quad(denom,0, len(TheoSpec)-1, (TheoSpec,ExperSpec,eTheo),limit =100))[0])
def ypendrySub(TheoSpec,ExperSpec):
    eTheo = ImaginaryEnergy(TheoSpec)
    #eExper = ImaginaryEnergy(ExperSpec)
    #specDif = SpecDifferenceSq(TheoSpec,ExperSpec)
    sumPendry = 0
    
    for n in range(len(TheoSpec)-1):
       sumPendry+=num(n,TheoSpec,ExperSpec,eTheo )/denom(n,TheoSpec,ExperSpec,eTheo)
        
    return sumPendry/4000
def Ypendry3(TheoSpec,ExperSpec):
    eTheo = ImaginaryEnergy(TheoSpec)
    sumPendry =0
    for wavelength in range(len(TheoSpec)-1):
       sumPendry+=(newY(wavelength, TheoSpec,eTheo) - newY(wavelength, ExperSpec,eTheo))**2 \
       / ((newY(wavelength, TheoSpec,eTheo)**2) + (newY(wavelength, ExperSpec,eTheo)**2))
    return sumPendry/4000
def file2Spectra(path):
    #Open and read file
    logfile =  open(path, 'r')
    logtest = logfile.read()
    logfile.close()
    # Find Frequencies with diceimal digits
    freqstri =re.findall('Frequencies\D*(\d+.\d+)\D*(\d+.\d+)\D*(\d+.\d+)',logtest) #looking for decimal numbers and spaces
    IRIntenstri =re.findall('IR Inten\D*(\d+.\d+)\D*(\d+.\d+)\D*(\d+.\d+)',logtest)
    IrDict =[]
    for freqTuple,intTuple in zip(freqstri,IRIntenstri):
        for n,p in zip(freqTuple,intTuple):
            IrDict.append( [float(n), float(p)])
    
    Irs = np.array(IrDict)
    #normalize
    Irs[:,1] = 100*Irs[:,1]/np.amax(Irs[:,1])
    return Irs
def gaussian_broadening(spectra, broaden, resolution=1):
 
    """ Performs gaussian broadening on IR spectrum
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """

    IR = np.zeros((int(4000/resolution) + 1))
    X = np.linspace(0,4000, int(4000/resolution)+1)
   # for f, i in zip(spectra[:,0], :  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
   # self.IR=np.vstack((X, IR)).T #tspec
                    
                    
    for line in spectra:
        freq = line[0]
        inten = line[1]
        IR += inten*np.exp(-0.5*((X-freq)/int(broaden))**2)
    return IR
fileList = glob.glob('Tri_A1*/Tri_A1*/input.log')
IR0Name =random.choice(fileList)
IR0 = file2Spectra(IR0Name)
IRBroad = gaussian_broadening(IR0,25,1)
Energy0 = ImaginaryEnergy(IRBroad)
IR1Name =random.choice(fileList)
IR1 = file2Spectra(IR1Name)
IRBroad1 = gaussian_broadening(IR1,25,1)
plt.plot(IRBroad, color= 'red')
plt.plot(IRBroad1, color= 'green')
Energy1 = ImaginaryEnergy(IRBroad1)
numJ =np.zeros(4000)
denomJ = np.zeros(4000)
for n in range(4000):
    numJ[n]= num(n,IRBroad, IRBroad1,Energy0 )
    denomJ[n] =num(n,IRBroad, IRBroad1,Energy0 )
#print(scipy.stats.wasserstein_distance(IRBroad,IRBroad/2))
#print(np.sum(numJ)-np.sum(denomJ))

print(ypendry(IRBroad, IRBroad1))
print(ypendrySub(IRBroad, IRBroad1))
print(Ypendry3(IRBroad, IRBroad))
timeB = time.perf_counter_ns()
#print((timeB-timeA)/10**9)
#print(ypendry(IRBroad,IRBroad+1))