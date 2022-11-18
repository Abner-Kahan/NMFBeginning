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





def nmfMatcher(OG_spectra,Calc_spectra):
    #pdb.set_trace() 
        #print(len(OG_spectra))
    #OG_spectra = np.transpose(OG_spectra)
    mindim = np.amin(OG_spectra.shape)
    errorTable = np.zeros((mindim,mindim))
    for n in range (mindim):
         for p in range(mindim):
                 errorTable[n,p] += np.sum(abs( OG_spectra[n] - Calc_spectra[p]))
    #print("hi \n", errorTable)
   # print(errorTable)
    matchTable=np.zeros((mindim,mindim))
    #print("errorTable \n \n",errorTable)
    for entry in range(mindim):
         Match = np.where(np.amin(errorTable) == errorTable)
         Match = list(Match) 
         
        
         x = Match[0]
         y = Match[1]
         matchTable[entry,0] =  x
         matchTable[entry,1] =  y
         #print(Match, errorTable[Match])
         errorTable[x]=10**12
        # errorTable[1,y]=10**12
         #print(errorTable)
         #print(errorTable)
    matchTable = matchTable.astype(int)
    return matchTable


fileList = glob.glob('Tri_A1*/Tri_A1*/input.log')

#print(isinstance(fileList[0],str))
#plt.vlines(file2Spectra(fileList[1])[:,0],0,file2Spectra(fileList[0])[:,1])


# In[52]:


def VertPlotParamaters():
    plt.xlabel("Wavenumber cm^-1")
    plt.ylabel("Intensity %")
    plt.gca().invert_yaxis() 
    plt.gca().invert_xaxis()
    plt.show()
    plt.clf()
    


# In[53]:

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
    return ( ( scipy.integrate.quad(num,0, len(TheoSpec)-1, (TheoSpec,ExperSpec, ImaginaryEnergy(TheoSpec)))[0]) / (
    scipy.integrate.quad(denom,0, len(TheoSpec)-1, (TheoSpec,ExperSpec,ImaginaryEnergy(TheoSpec))))[0])
def nmf2TesterMixB():
    #pdb.set_trace()
    fig, axs = plt.subplots(3,2)
    fraction1 = random.random()
    fraction2 = random.random()
    print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}.')
    
    #Creating Two Spectra
    IR0Name =random.choice(fileList)
    print(IR0Name)
    IR0 = file2Spectra(IR0Name)
    IR0 = gaussian_broadening(IR0,25,1)
    axs[0,0].plot(IR0/np.max(IR0), color= 'blue')
    axs[0,0].set_title("Glycan 1", x=.5, y= .75)


    
    IR1Name =random.choice(fileList)
    print(IR1Name)
    IR1 = file2Spectra(IR1Name)
    IR1 = gaussian_broadening(IR1,25,1)
    axs[0,1].plot(IR1/np.max(IR1), color= 'red')
    axs[0,1].set_title("Glycan 2", x=.5, y= .75)
   # print(IR1.shape)
    print("Correlation between Functions", ypendry(IR0,IR1))
    IRF = np.zeros((2,4001))
    IRF[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRF[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    axs[1,0].plot(IRF[0,:]/np.max(IRF[0,:]), color='purple')
    axs[1,0].set_title("Mixture 1", x=.5, y= .75)
    axs[1,1].plot(IRF[1,:]/np.max(IRF[1,:]), color = 'purple')
    axs[1,1].set_title("Mixture 2", x=.5, y= .75)
    
    IRF= np.transpose(IRF)
    
    model = NMF(n_components=2, max_iter=4000, tol= 1*10**-10, solver= 'mu', init ='nndsvda', beta_loss= 'kullback-leibler' )
    W = model.fit_transform(IRF)
    H = model.components_
    Hnorm = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    print("H", H)
    print("H_Norm", Hnorm)
    product = np.matmul(W,H)
    matchTable = nmfMatcher (IRF, W)
    IrOrgs = [IR0,IR1]
    #axs[4].plot(product[:,matchTable[0,1]])
    

    
    
    scale1  = np.mean(IR0)/np.mean(W[:,0])
    scale2  = np.mean(IR1)/np.mean(W[:,0])
    #difference between first input and first NMF
    difA = np.sum(abs(IR0 - scale1 * W[:,0] ))
    difB = np.sum(abs(IR1 - scale2 * W[:,0] ))
    if difA < difB:

        axs[2,0].plot(W[:,0]/np.max(W[:,0]), color= 'blue')
        
        axs[2,1].plot(W[:,1]/np.max(W[:,1]), color = 'red')
    else:
        axs[2,0].plot(W[:,1]/np.max(W[:,1]), color= 'blue' )
        
        axs[2,1].plot(W[:,0]/np.max(W[:,0]), color = 'red')  
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
    print(matchTable)
  #  print(ypendry(IrOrgs[matchTable[0,0]], product[:,matchTable[0,1]]))
   # print(ypendry(IrOrgs[matchTable[1,0]], product[:,matchTable[1,1]]))
   # print("Matrix product size", product.shape)
    ##print("Variable matrix", H)
   # print("Product Sizes", W.shape, IRF.shape)
    #print(nmfMatcher (IRF, product))
 
nmf2TesterMixB()






