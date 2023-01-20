#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import chirp, find_peaks, peak_widths
import scipy

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
def percentError(Original, NMF):
    return 100.0*(NMF-Original)/Original
def getFrac(aray):
    aray[aray == 0.0] = 10**-16
    k = (aray[1,0]/aray[1,1] - 1.0 )/    \
    ( aray[0,1] * aray[1,0] / ( aray[0,0] * aray[1,1] ) - 1.0)
    j = ( aray[1,0]/aray[1,1] -1) / ( aray[1,0]/aray[1,1]- aray[0,0]/aray[0,1] )
    return k,j

def ImaginaryEnergy(spectra):
    peaks, _ = find_peaks(spectra)
    results_half = peak_widths(spectra, peaks, rel_height=0.5)
    ImaginaryEnergy = np.average (results_half[0])/2.0
    return ImaginaryEnergy


def newY(wavelength, spectra,V_oi):
    if spectra[wavelength] == 0:
        spectra[wavelength] = 10**-16
    #L = scipy.misc.derivative (lambda x:(spectra[int(x)]), wavelength)/spectra[int(wavelength)]
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

#nndsvd


# In[3]:


def IrPlotter(item,title,leg = [], multiple = False):
    if not(multiple):
        plt.plot(np.linspace(0,4000,len(item)),item,markersize=.1)
    else:
        for n in item:
            plt.plot(np.linspace(0,4000,len(n)),n,markersize=.1)
    if len (leg) > 0:
        plt.legend(leg)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")    
    plt.show()
    plt.clf()

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

#print(isinstance(fileList[0],str))
#plt.vlines(file2Spectra(fileList[1])[:,0],0,file2Spectra(fileList[0])[:,1])




# In[53]:



def nmf2TesterMixB():
    #pdb.set_trace()
    fig, axs = plt.subplots(3,2)
    fraction1 = random.random()
    fraction2 = random.random()
    print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}.')
    IR0 = addIRS(8,1000)
    IR1 = addIRS(8,1000)

    axs[0,0].plot(IR0/np.max(IR0), color= 'goldenrod')
    axs[0,0].set_title("RS 1", x=.5, y= .75)


    axs[0,1].plot(IR1/np.max(IR1), color= 'dodgerblue')
    axs[0,1].set_title("RS 2", x=.5, y= .75)
   # print(IR1.shape)
    #print("Correlation between Functions", ypendry(IR0,IR1))
    IRF = np.zeros((2,1000))
    IRFPure = np.zeros((4,1000))
    IRF[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRF[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    axs[1,0].plot(IRF[0,:]/np.max(IRF[0,:]), color='olivedrab')
    axs[1,0].set_title("Mixture 1", x=.5, y= .75)
    axs[1,1].plot(IRF[1,:]/np.max(IRF[1,:]), color = 'olivedrab')
    axs[1,1].set_title("Mixture 2", x=.5, y= .75)
    
    IRFPure[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRFPure[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    IRFPure[2,:] = IR0
    IRFPure[3,:] = IR1
    
    
    
    IRF= np.transpose(IRF)
    IRFPure = np.transpose(IRFPure)
    model = NMF(n_components=2, max_iter=1000, tol= 1*10**-10, solver= 'mu', init ='nndsvda', beta_loss= 'kullback-leibler' )
    W = model.fit_transform(IRF)
    H = model.components_
    
    WPure = model.fit_transform(IRFPure)
    HPure = model.components_
    Hnorm = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    #print("H", H)
   # print("H_Norm", Hnorm)
    product = np.matmul(W,H)
    IrOrgs = [IR0,IR1]
    #axs[4].plot(product[:,matchTable[0,1]])
    
    
    
    #scale1  = np.mean(IR0)/np.mean(WPure[:,0])
   # scale2  = np.mean(IR0)/np.mean(WPure[:,1])
    #difference between first input and first NMF
   # RealdifA = np.sum(abs(IR0 - scale1 * WPure[:,0] ))
    #RealdifB = np.sum(abs(IR0 - scale2 * WPure[:,1] ))
    ypendryPure = Ypendry3(IR0,WPure[:,0])
    ypendryPureFlip = Ypendry3(IR0,WPure[:,1])
    if ypendryPureFlip < ypendryPure: 
        HPure[[0, 1]] = HPure[[1, 0]] 
        WPure[:, [1, 0]] = WPure[:, [0, 1]]
        ypendryPure = ypendryPureFlip
   # print(HPure)

    ypendry = Ypendry3(IR0,W[:,0])
    ypendry2 = Ypendry3(IR0,W[:,1])
    ypendryFlip = Ypendry3(IR0,W[:,1])
    ypendryFlip2 = Ypendry3(IR0,W[:,1])
    if ypendryFlip + ypendryFlip2 < ypendry + ypendry2: 
        H[[0, 1]] = H[[1, 0]] 
        W[:, [1, 0]] = W[:, [0, 1]]
        ypendry = ypendryFlip
        print("Flip \n\n\n\n\n")
        
        
        
    axs[2,0].plot(W[:,0]/np.max(W[:,0]), color= 'goldenrod')
        
    axs[2,1].plot(W[:,1]/np.max(W[:,1]), color = 'dodgerblue')
   
        
        
        
        
        
    axs[2,0].set_title("NMF RS 1", x=.5, y= .75)
    axs[2,1].set_title("NMF RS 2", x=.5, y= .75)
    #IrPlotter([IrOrgs[matchTable[0,0]], product[:,matchTable[0,1]]],"First Calculated Spectra", [IR0Name[9:-10], 'NMF Generated IR'],True)
    #IrPlotter([IrOrgs[matchTable[1,0]] ,product[:,matchTable[1,1]]],"Second Calculated Spectra", [IR1Name[9:-10], 'NMF Generated IR'], True)      
   # IrPlotter([IROrgs[matchTable[0][0]],product
               #[W[matchTable[0][1],:]]], "First Matched Spectra", True)
    axs[2,0].set(xlabel='cm$^-1$')
    axs[2,1].set(xlabel='cm$^-1$')
    axs[0,0].set(xticks = [])
    axs[1,0].set(xticks = [])
    axs[0,1].set(xticks = [])
    axs[1,1].set(xticks = [])
    plt.show()
    plt.clf()
    fig2, axs2 = plt.subplots(2,2)
    axs2[0,0].plot(IR0/np.max(IR0))
    axs2[0,0].set_title("IR0")
    axs2[0,1].plot(IR1/np.max(IR1))
    axs2[0,1].set_title("IR1")
    axs2[1,0].plot(WPure[:,0]/np.max(WPure[:,0]))
    axs2[1,0].set_title("NMF 0")
   
    axs2[1,1].plot(WPure[:,1]/np.max(WPure[:,1]))
    axs2[1,1].set_title("NMF 1")
   # plt.show()
   # plt.plot(W[:,0])
   # plt.show()
    product = np.matmul(W,H)
    print(product.shape)
    print(W.shape)
    print(WPure.shape)
    PureProduct = np.matmul(WPure,HPure)
    print("Impure \n")
    print (f"The expected fractions from the mixture are {getFrac(H)[0]:.6} and {getFrac(H)[1]:.6}")    
    print (f"Percent error is {percentError(fraction1,getFrac(H)[0])}  \
           and {percentError(fraction2,getFrac(H)[1])}")
    print (f"The expected fractions from the mixture are {1-getFrac(H)[0]:.6} and {1-getFrac(H)[1]:.6}")    
    print (f"Percent error is {percentError(1-fraction1,1-getFrac(H)[0])}  \
               and {percentError(1-fraction2,1-getFrac(H)[1])}")   
           
    print(f"Ypendry error is {ypendry}")
    print(f"Ypendry error is {Ypendry3(IR1,W[:,1])}")    
    print("Pure \n")
    print (f"The expected fractions from the mixture are {getFrac(HPure)[0]:.6} and {getFrac(HPure)[1]:.6}")    
    print (f"Percent error is {percentError(fraction1,getFrac(HPure)[0])}  \
           and {percentError(fraction2,getFrac(HPure)[1])}")
    print (f"Percent error is {percentError(1-fraction1,1-getFrac(HPure)[0])}  \
           and {percentError(1-fraction2,1-getFrac(HPure)[1])}")      
    print(f"Ypendry error is {ypendryPure}")
    print(f"Ypendry error is {Ypendry3(IR1,WPure[:,1])}")    
    print(f"Ypendry error of two spectra is  {Ypendry3(IR0,IR1)}")    
  #  print(ypendry(IrOrgs[matchTable[0,0]], product[:,matchTable[0,1]]))
   # print(ypendry(IrOrgs[matchTable[1,0]], product[:,matchTable[1,1]]))
   # print("Matrix product size", product.shape)
    ##print("Variable matrix", H)
   # print("Product Sizes", W.shape, IRF.shape)
    #print(nmfMatcher (IRF, product))

    
   
 
nmf2TesterMixB()

def addIRS2(peakNum,PointNum,xloc):
    IR = np.zeros((PointNum,))
    X =  np.linspace(0, PointNum, PointNum)
    if type(xloc) ==int:
        xloc =[xloc]
    peakHeigh = [ .6 for i in range(peakNum) ] 
    thickness = [ PointNum * .02 for i in range(peakNum) ] 

    for f, i, b in zip(xloc, peakHeigh, thickness):  
        IR += i*np.exp(-0.5*((X-f)/int(b))**2)
    return IR

plt.clf()
# =============================================================================
# fig3, axs3 = plt.subplots(3)
# spec1 =addIRS2(1,1000,200)
# 
# spec2 =addIRS2(1,1000,240)
# MixedSpec =   addIRS2(2,1000,[200,240]) 
# axs3[0].plot(spec1/np.max(MixedSpec), color='red')
# axs3[1].plot(spec2/np.max(MixedSpec), color='blue') 
# axs3[0].set(ylim = [0,1])
# axs3[1].set(ylim = [0,1])
# axs3[0].set(xticks = [])
# axs3[1].set(xticks = [])
# axs3[2].plot(MixedSpec/np.max(MixedSpec), color='darkviolet')
# axs3[2].set(xlabel='cm$^-1$')
# axs3[0].set(title='Spectra 1')
# axs3[1].set(title='Spectra 2')
# axs3[2].set(title='Mixed IR Spectra')
# plt.show(fig3)
# =============================================================================
