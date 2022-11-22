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
    
    #Creating Two Spectra
    IR0Name =random.choice(fileList)
  #  print(IR0Name)
    IR0 = file2Spectra(IR0Name)
    IR0 = gaussian_broadening(IR0,25,1)
    axs[0,0].plot(IR0/np.max(IR0), color= 'orangered')
    axs[0,0].set_title("Glycan 1", x=.5, y= .75)


    
    IR1Name =random.choice(fileList)
    print(IR1Name)
    IR1 = file2Spectra(IR1Name)
    IR1 = gaussian_broadening(IR1,25,1)
    axs[0,1].plot(IR1/np.max(IR1), color= 'deepskyblue')
    axs[0,1].set_title("Glycan 2", x=.5, y= .75)
   # print(IR1.shape)
    #print("Correlation between Functions", ypendry(IR0,IR1))
    IRF = np.zeros((2,4001))
    IRFPure = np.zeros((4,4001))
    IRF[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRF[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    axs[1,0].plot(IRF[0,:]/np.max(IRF[0,:]), color='darkorchid')
    axs[1,0].set_title("Mixture 1", x=.5, y= .75)
    axs[1,1].plot(IRF[1,:]/np.max(IRF[1,:]), color = 'darkorchid')
    axs[1,1].set_title("Mixture 2", x=.5, y= .75)
    
    IRFPure[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRFPure[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    IRFPure[2,:] = IR0
    IRFPure[3,:] = IR1
    
    
    
    IRF= np.transpose(IRF)
    IRFPure = np.transpose(IRFPure)
    model = NMF(n_components=2, max_iter=4000, tol= 1*10**-10, solver= 'mu', init ='nndsvda', beta_loss= 'kullback-leibler' )
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
    

    
    
    scale1  = np.mean(IR0)/np.mean(WPure[:,0])
    scale2  = np.mean(IR0)/np.mean(WPure[:,1])
    #difference between first input and first NMF
    RealdifA = np.sum(abs(IR0 - scale1 * WPure[:,0] ))
    RealdifB = np.sum(abs(IR0 - scale2 * WPure[:,1] ))
    if RealdifA > RealdifB: 
        HPure[[0, 1]] = HPure[[1, 0]] 
        WPure[:, [1, 0]] = WPure[:, [0, 1]]
   # print(HPure)
    scale1  = np.mean(IR0)/np.mean(W[:,0])
    scale2  = np.mean(IR0)/np.mean(W[:,1])
 
    
    #difference between first input and first NMF
    difA = np.sum(abs(IR0 - scale1 * W[:,0] ))
    difB = np.sum(abs(IR0 - scale2 * W[:,1] ))
    
    
    if difA < difB:
        H[[0, 1]] = H[[1, 0]] 
        W[:, [1, 0]] = W[:, [0, 1]]
        
        
        
    axs[2,0].plot(W[:,0]/np.max(W[:,0]), color= 'orangered')
        
    axs[2,1].plot(W[:,1]/np.max(W[:,1]), color = 'deepskyblue')
   
        
        
        
        
        
    axs[2,0].set_title("NMF Glycan 1", x=.5, y= .75)
    axs[2,1].set_title("NMF Glycan 2", x=.5, y= .75)
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
    plt.plot(IR0)
    plt.show()
    plt.plot(W[:,0])
    plt.show()
    product = np.matmul(W,H)
    print(product.shape)
    print(W.shape)
    print(WPure.shape)
    PureProduct = np.matmul(WPure,HPure)
    print("Impure \n")
    print (f"The expected fractions from the mixture are {getFrac(H)[0]:.6} and {getFrac(H)[1]:.6}")    
    print (f"Percent error is {percentError(fraction1,getFrac(H)[0])}  \
           and {percentError(fraction2,getFrac(H)[1])}")
    print(f"Ypendry error is {Ypendry3(IR0,W[:,0])}")
    print(f"Ypendry error is {Ypendry3(IR1,W[:,1])}")    
    print("Pure \n")
    print (f"The expected fractions from the mixture are {getFrac(HPure)[0]:.6} and {getFrac(HPure)[1]:.6}")    
    print (f"Percent error is {percentError(fraction1,getFrac(HPure)[0])}  \
           and {percentError(fraction2,getFrac(HPure)[1])}")
    print(f"Ypendry error is {Ypendry3(IR0,WPure[:,0])}")
    print(f"Ypendry error is {Ypendry3(IR1,WPure[:,1])}")    
    print(f"Ypendry error of two spectra is  {Ypendry3(IR0,IR1)}")    
  #  print(ypendry(IrOrgs[matchTable[0,0]], product[:,matchTable[0,1]]))
   # print(ypendry(IrOrgs[matchTable[1,0]], product[:,matchTable[1,1]]))
   # print("Matrix product size", product.shape)
    ##print("Variable matrix", H)
   # print("Product Sizes", W.shape, IRF.shape)
    #print(nmfMatcher (IRF, product))
 
nmf2TesterMixB()






