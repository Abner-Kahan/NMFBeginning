#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:27:04 2022

@author: abnerkahan
"""

#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF

#np.set_printoptions(precision=3)



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
    print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}.')
    fig, axs = plt.subplots(3,2)

    #unscaled
    IRF = np.zeros((2,numPoints))
    IR0 = addIRS(8,numPoints)
    axs[0,0].plot(IR0/np.max(IR0), color= 'goldenrod')
    axs[0,0].set_title("RS 1", x=.5, y= .75)
    
  
    #IrPlotter(IR0,"FirstPlot")
    IR1 = addIRS(8,numPoints)
    axs[0,1].plot(IR1/np.max(IR1), color= 'dodgerblue')
    axs[0,1].set_title("RS 2", x=.5, y= .75)

    
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((4,numPoints))
    IrMix[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMix[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    IrMix[2,:] = IR0
    IrMix[3,:] = IR1  
    #normalized
    
    #IR Mix less data
    IrMixLD = np.zeros((2,numPoints))
    IrMixLD[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMixLD[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    

    #IrMix[3,:] = IR0*(1-fraction2)  + IR1*fraction1
    axs[1,0].plot(IrMix[0,:]/np.max(IrMix[0,:]), color= 'olivedrab')
    axs[1,0].set_title("Mixture 1", x=.5, y= .75)
    axs[1,1].plot(IrMix[1,:]/np.max(IrMix[1,:]), color= 'olivedrab')
    axs[1,1].set_title("Mixture 2", x=.5, y= .75)
    

   # customH = 
    IrMixLD = np.transpose(IrMixLD)
    IrMix= np.transpose(IrMix)
   # model = NMF(n_components=2, init='nndsvd',  max_iter=15000, tol= 1*10**-7)
    model = NMF(n_components=2, max_iter=1000, tol= 1*10**-10, solver= 'mu', init ='nndsvda', beta_loss= 'kullback-leibler' )
    

    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    H = model.components_
    W_ld = model.fit_transform(IrMixLD)
    H_ld = model.components_
    

    print(H)
   # print("bob")
    H0 = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
   # print(H0)
    #print("bob2")
    #print(H_ld) 
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
    
   # print (difA, difB)
   # print (RealdifA, RealdifB)
    
    if difA < difB:

        axs[2,0].plot(np.linspace(0,1000,10000), W_ld[:,0]/np.max(W_ld[:,0]), color= 'goldenrod')
        axs[2,1].plot(np.linspace(0,1000,10000), W_ld[:,1]/np.max(W_ld[:,1]), color= 'dodgerblue')
        
       
        # print ("Regular")
        # print (H[0,0]/(H[1,0]+H[0,0]))
        # print(H[0,1]/(H[1,1]+H[0,1]))   
        
        # print ("Norm")
        # print (H[0,0]/np.sum(IR0)/(H[0,1]/np.sum(IR1)+H[0,0]/np.sum(IR0)))
        # print(H[1,0]/np.sum(IR0)/(H[1,1]/np.sum(IR1)+H[1,0]/np.sum(IR0)))
        

        
 
    else:
        print("flip\n\n")
        axs[2,0].plot(np.linspace(0,1000,10000), W_ld[:,1]/np.max(W_ld[:,1]), color= 'goldenrod')
        axs[2,1].plot(np.linspace(0,1000,10000), W_ld[:,0]/np.max(W_ld[:,0]), color= 'dodgerblue')
        H_ld[[0, 1]] = H_ld[[1, 0]]  
        
    print (f"The expected fractions from the mixture are {getFrac(H_ld)[0]:.6} and {getFrac(H_ld)[1]:.6}")    
    print (f"Percent error is {(getFrac(H_ld)[0]-fraction1)/fraction1 *100}  \
           and {(getFrac(H_ld)[1]-fraction2)/fraction2 *100}")
           
    if RealdifA > RealdifB: 
        H[[0, 1]] = H[[1, 0]] 
        print ("Real Flip \n\n")
    print (f"The expected fractions are {H[0,0]/np.max(H[0]):.6} and {H[0,1]/np.max(H[0]):.6}")
        
    print ((f"The percent error is  {(H[0,0]/np.max(H[0]) -fraction1 )/fraction1*100:.6} and  \
                {(H[0,1]/np.max(H[0]) -fraction2)/fraction2*100:.6}"))
   
    
    
    
    
  #  print (f"The expected fractions are {H_ld[0,0]/np.sum(H_ld[0]):.6} and {H_ld[1,0]/np.sum(H_ld[1]):.6}")
   # print (f"The expected fractions are {H_ld[0,0]/(H_ld[0,0] + H_ld[1,0] ):.6} and {H_ld[0,1]/(H_ld[1,1] + H_ld[0,1] ):.6}")
          
        # print ("Regular")        
        # print(H[1,0]/np.sum(IR0)/(H[1,1]/np.sum(IR1)+H[1,0]/np.sum(IR0)))
        # print (H[0,0]/np.sum(IR0)/(H[0,1]/np.sum(IR1)+H[0,0]/np.sum(IR0)))
        # print ("Norm")           
        # print(H[1,0]/(H[1,1]+H[1,0]))
        # print (H[0,0]/(H[0,1]+H[0,0]))
        
        # print ("Norm NMF")

        

    
    axs[2,0].set_title("NMF RS 1 ", x=.5, y= .75)
    axs[2,1].set_title("NMF RS 2", x=.5, y= .75)
    
    
    
    
    axs[0,0].set(xticks = [])
    axs[1,0].set(xticks = [])
    axs[0,1].set(xticks = [])
    axs[1,1].set(xticks = [])
    axs[2,0].set(xlabel='cm$^-1$')
    axs[2,1].set(xlabel='cm$^-1$')
    
    
    

    plt.show(fig)
    #plt.show(figUS)
    return H
    #axs[2,0].set(xticks = [np.linspace(0, numPoints, 11).all()])
   # axs[2,1].set(xticks = [np.linspace(0, numPoints, 11).all()])
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
    #print(model.fit(IrMix))
# =============================================================================
    # for entri in nmfMatcher(IRF,W):
        
    #      plt.plot(np.linspace(0,1000,numPoints),IRF[:,entri[0][0]],color="red")
    #      if H[0,0]>.01:
    #          print("The fraction of the first is", H[0,2])
    #          print("The fraction of the second is", H[1,2])

             
    #      else:
    #          print("The fraction of the first is", H[1,2])
    #          print("The fraction of the second is", H[0,2])
    #      plt.plot(np.linspace(0,1000,numPoints),(W[:,entri[1][0]]*(max(HO[entri[1][0]]))))
    #     # print("full", (max(HO[entri[0][0]])))
        


    #      plt.gca().invert_yaxis()
    #      plt.legend(["Original", "Calculated"])
    #      plt.title(str(entri[0][0])+ " Both Spectra")
    #      plt.xlabel("cm^-1")
    #      plt.show()
    #      plt.clf()
# =============================================================================


#nmfTesterMix([1,.5,.25,0,.8,.3],10000)
# =============================================================================
H = nmf2TesterMix(10000)


