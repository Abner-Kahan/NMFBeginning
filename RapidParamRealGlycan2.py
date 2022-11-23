#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 17:25:09 2022

@author: abnerkahansmack
"""
import glob
import re
import sys


import time
import warnings
#warnings.filterwarnings("ignore", message="/home/abnerkahan/anaconda3/lib/python3.9/site-packages/sklearn/decomposition/_nmf.py:1637: ConvergenceWarning:e")
#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF
from scipy.signal import chirp, find_peaks, peak_widths
import scipy

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
    return ( ( scipy.integrate.quad(num,0, len(TheoSpec)-1, (TheoSpec,ExperSpec, ImaginaryEnergy(TheoSpec)),limit =100)[0]) / (
    scipy.integrate.quad(denom,0, len(TheoSpec)-1, (TheoSpec,ExperSpec,ImaginaryEnergy(TheoSpec)),limit =100))[0])
def percentError(Original, NMF):
    if Original ==0:
        print("hell nah")
    return abs(100.0*(NMF-Original)/Original)
def getFrac(aray):
    aray[aray == 0.0] = 10**-16
    k = (aray[1,0]/aray[1,1] - 1.0 )/    \
    ( aray[0,1] * aray[1,0] / ( aray[0,0] * aray[1,1] ) - 1.0)
    j = ( aray[1,0]/aray[1,1] -1) / ( aray[1,0]/aray[1,1]- aray[0,0]/aray[0,1] )
    
    #print(k,j)
    return k,j
def newY(wavelength, spectra,V_oi):
    if spectra[wavelength] == 0:
        spectra[wavelength] = 10**-16
    L = (spectra[wavelength+1]-spectra[wavelength])/spectra[wavelength]
    if L == 0 :
        L = 10**-16

    return (L**-1)/(L**-2 + V_oi**2)
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
paramlist = ['random_cd_frobenius' ,
               'nndsvd_cd_frobenius' ,
               'nndsvda_cd_frobenius' ,
               'nndsvdar_cd_frobenius' , 
               'random_mu_frobenius' ,
               'nndsvd_mu_frobenius' ,
               'nndsvda_mu_frobenius' ,
               'nndsvdar_mu_frobenius' ,
               'random_mu_kullback-leibler' ,
               'nndsvd_mu_kullback-leibler' ,
               'nndsvda_mu_kullback-leibler' ,
               'nndsvdar_mu_kullback-leibler'  ]
def nmf2TesterMix(numPoints):

    fraction1 = random.random()
    fraction2 = random.random()
    IR0Name =random.choice(fileList)
   # print(IR0Name)
    IR0 = file2Spectra(IR0Name)
    IR0 = gaussian_broadening(IR0,25,1)

    IR1Name =random.choice(fileList)
   # print(IR1Name)
    IR1 = file2Spectra(IR1Name)
    IR1 = gaussian_broadening(IR1,25,1)

    
    # IrMix = np.zeros((4,numPoints))
    # IrMix[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    # IrMix[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    # IrMix[2,:] = IR0
    # IrMix[3,:] = IR1  
    
    #plt.plot(IR0)
   # plt.plot(IR1)
   # plt.show()
    #plt.clf()
    IrMixLD = np.zeros((2,numPoints))
    IrMixLD[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMixLD[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    IrOrg = IrMixLD
    IrMixLD = np.transpose(IrMixLD)
  #  IrMix= np.transpose(IrMix)
    ErrorDict ={}
    for param in paramlist:
        model = NMF(n_components=2, max_iter=1000, tol= 1*10**-10, \
        init =param[:param.find('_')], solver = param[param.find('_')+1:param.rfind('_')],  
        beta_loss = param[param.rfind('_')+1:])
            
        W_ld = model.fit_transform(IrMixLD)
        H_ld = model.components_
        Errt = model.reconstruction_err_       
        ypendry = Ypendry3(IR0,W_ld[:,0])
        ypendry2 = Ypendry3(IR0,W_ld[:,1])
        ypendryFlip = Ypendry3(IR0,W_ld[:,1])
        ypendryFlip2 = Ypendry3(IR0,W_ld[:,1])
        if ypendryFlip + ypendryFlip2 < ypendry + ypendry2:
           # print("flip")
           # plt.plot(IR0)
           # plt.plot(W_ld[:,1])
           # plt.show()
           # plt.clf()
            H_ld[[0, 1]] = H_ld[[1, 0]]
            W_ld[:, [1, 0]] = W_ld[:, [0, 1]]
            ypendry = ypendryFlip
            ypendry2 = ypendryFlip2
        NewFrac1 =  getFrac(H_ld)[0]
        NewFrac2 =  getFrac(H_ld)[1]
    
        PendryR = (ypendry+ypendry2)/2
     #   RealFrac1 = H[0,0]/np.max(H[0])
     #   RealFrac2 = H[0,1]/np.max(H[0])
       # if RealdifB < RealdifA:
       #     print("flip")
       #     RealFrac1 = 1 -RealFrac1
        #    RealFrac2 = 1- RealFrac2
       # wassertein1 = scipy.stats.wasserstein_distance(IrMixLD[:,0],product[:,0])
        #wassertein2 = scipy.stats.wasserstein_distance(IrMixLD[:,1],product[:,1])
    
        FractionError =    (percentError(fraction1,NewFrac1 ) + percentError(1.0 - fraction1,1.0- NewFrac1 ) +  \
            percentError(fraction2, NewFrac2) + percentError(1.0- fraction2, 1.0 - NewFrac2)  )/4.0                   
        # if Errt> 50:
        #     plt.plot(IR0)
        #     plt.plot(W_ld[:,0])
        #     plt.show()
        #     plt.clf()
        ErrorDict[str(param)] = [FractionError,Errt, PendryR]
    return ErrorDict

            #5,7,9,11
            #13,15,17,19
            
            
            
# [Fraction1, 1-Fraction, Fraction2, 1-Fraction2, 
# Frac1 -NMF, % error, 1-Frac1 NMF, %error, 
#Frac2 -NMF, % error, 1-Frac2 NMF, %error,
# Frac1 -NMF Real, % error, 1-Frac1 NMF Real, %error,
# Frac2 -NMF Real, % error, 1-Frac2 NMF Real, %error  ]    
timeA = time.perf_counter_ns() 

iters = 250
ResultsTable = pd.DataFrame(columns = paramlist)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for n in range(iters):
        ResultsTable = pd.concat([ResultsTable, pd.DataFrame(nmf2TesterMix(4001))])
print ("% Error")
print(ResultsTable.loc[0,:].mean())
print ("Fro Error")
print(ResultsTable.loc[1,:].mean())
print ("Pendry Error")
print(ResultsTable.loc[2,:].mean())
    


#print(ResultsTable)
#BigTable2 = ResultsTable[~np.isnan(ResultsTable).any(axis=1), :]
#print ("Shape", BigTable2.shape)
#print(ResultsTable)

#print(BigTable2)
#print(np.mean(abs(BigTable2[:,0])))
timeB = time.perf_counter_ns()
print((timeB-timeA)/10**9)
sys.stdout.write('\a')
sys.stdout.flush()
#print(np.mean(abs(BigTable2[:,1])))
#print(np.mean(BigTable2[:,2]))

#print(np.mean(BigTable2[:,3]))
#np.save('BigIterFractions.npy',ResultsTable)

    
