#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 17:25:09 2022

@author: abnerkahansmack
"""

#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF


def percentError(Original, NMF):
    return 100*(NMF-Original)/Original
def getFrac(aray):
    k = (aray[1,0]/aray[1,1] - 1 )/    \
    ( aray[0,1] * aray[1,0] / ( aray[0,0] * aray[1,1] ) - 1)
    j = (1- aray[1,0]/aray[1,1]) / (aray[0,0]/aray[0,1] - aray[1,0]/aray[1,1])
    return k,j
def addIRS(peakNum,PointNum):
    IR = np.zeros((PointNum,))
    X =  np.linspace(0, PointNum, PointNum)

    xloc = [ random.randrange(0, PointNum) for i in range(peakNum) ] 
    peakHeigh = [ random.random() for i in range(peakNum) ] 
    thickness = [ PointNum * random.choice([.007,.008,.009,.01,.011,.012,.013,.014,.042]) for i in range(peakNum) ] 

    for f, i, b in zip(xloc, peakHeigh, thickness):  
        IR += i*np.exp(-0.5*((X-f)/int(b))**2)

    plt.clf
    return IR 

def nmf2TesterMix(numPoints):
    fraction1 = random.random()
    fraction2 = random.random()

    IR0 = addIRS(8,numPoints)
    IR1 = addIRS(8,numPoints)
    
    IrMix = np.zeros((4,numPoints))
    IrMix[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMix[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    IrMix[2,:] = IR0
    IrMix[3,:] = IR1  
    
    
    IrMixLD = np.zeros((2,numPoints))
    IrMixLD[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMixLD[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    
    IrMixLD = np.transpose(IrMixLD)
    IrMix= np.transpose(IrMix)
    
    model = NMF(n_components=2, max_iter=800, tol= 1*10**-10, solver= 'mu', init ='nndsvda', beta_loss= 'kullback-leibler' )
    W = model.fit_transform(IrMix)
    H = model.components_
    W_ld = model.fit_transform(IrMixLD)
    H_ld = model.components_
    
    
    scale1  = np.mean(IR0)/np.mean(W_ld[:,0])
    scale2  = np.mean(IR0)/np.mean(W_ld[:,1])
    #difference between first input and first NMF
    difA = np.sum(abs(IR0 - scale1 * W_ld[:,0] ))
    difB = np.sum(abs(IR0 - scale2 * W_ld[:,1] ))
    
    Realscale1  = np.mean(IR0)/np.mean(W[:,0])
    Realscale2  = np.mean(IR0)/np.mean(W[:,1])
    #difference between first input and first NMF
    RealdifA = np.sum(abs(IR0 - scale1 * W[:,0] ))
    RealdifB = np.sum(abs(IR0 - scale2 * W[:,1] ))
    
    if difB < difA:
        H_ld[[0, 1]] = H_ld[[1, 0]]
        

    NewFrac1 = getFrac(H_ld)[0]
    NewFrac2 = getFrac(H_ld)[1]
    
    RealFrac1 = H[0,0]/np.max(H[0])
    RealFrac2 = H[0,1]/np.max(H[0])
    if RealdifB < RealdifA:
        print("flip")
        RealFrac1 = 1 -RealFrac1
        RealFrac2 = 1- RealFrac2
        
                              
    
    return fraction1, fraction2, 1-fraction1, 1-fraction2,  \
        NewFrac1, percentError(fraction1,NewFrac1 ),1-NewFrac1, percentError(1 - fraction1,1- NewFrac1 ),  \
        NewFrac2, percentError(fraction2, NewFrac2), 1- NewFrac2, percentError(1- fraction2, 1 - NewFrac2), \
        RealFrac1, percentError(fraction1,RealFrac1 ),1-RealFrac1, percentError(1 - fraction1,1- RealFrac1 ),  \
        RealFrac2, percentError(fraction2, RealFrac2), 1- RealFrac2, percentError(1- fraction2, 1 - RealFrac2)

            #5,7,9,11
            #13,15,17,19
            
            
            
# [Fraction1, 1-Fraction, Fraction2, 1-Fraction2, 
# Frac1 -NMF, % error, 1-Frac1 NMF, %error, 
#Frac2 -NMF, % error, 1-Frac2 NMF, %error,
# Frac1 -NMF Real, % error, 1-Frac1 NMF Real, %error,
# Frac2 -NMF Real, % error, 1-Frac2 NMF Real, %error  ]    
    
iters = 10
ResultsTable = np.zeros((iters, 20))
for n in range(iters):
    ResultsTable[n] =  nmf2TesterMix(10000)
    
np.save('BigIterFractions.npy',ResultsTable)

    