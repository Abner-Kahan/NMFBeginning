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

def file2Spectra(path):
    logfile =  open(path, 'r')
    logtest = logfile.read()
    logfile.close()
    freqstri =re.findall('Frequencies\D*(\d+.\d+)\D*(\d+.\d+)\D*(\d+.\d+)',logtest) #looking for decimal numbers and spaces
    IRIntenstri =re.findall('IR Inten\D*(\d+.\d+)\D*(\d+.\d+)\D*(\d+.\d+)',logtest)
    IrDict =[]
    for freqTuple,intTuple in zip(freqstri,IRIntenstri):
        for n,p in zip(freqTuple,intTuple):
            IrDict.append( [float(n), float(p)])
    
    Irs = np.array(IrDict)
    Irs[:,1] = 100*Irs[:,1]/np.amax(Irs[:,1])
    return Irs
# =============================================================================
# for Plotting
#     plt.vlines(Irs[:,0],0,Irs[:,1])
#     plt.xlabel("Wavenumber cm^-1")
#     plt.ylabel("Intensity %")
#     plt.gca().invert_yaxis() 
#     plt.gca().invert_xaxis() 
# =============================================================================




def nmfMatcher(OG_spectra,Calc_spectra):
 
        #print(len(OG_spectra))
    OG_spectra = np.transpose(OG_spectra)
    errorTable = np.zeros((OG_spectra.shape[1], Calc_spectra.shape[1]))
    for n in range (OG_spectra.shape[0]):
         for p in range(OG_spectra.shape[1]):
             for q in range(Calc_spectra.shape[1]):
                 errorTable[p,q] += abs( OG_spectra[n,p] - Calc_spectra[n,q])
    print("hi \n", errorTable)
    matchTable=[]
    #print("errorTable \n \n",errorTable)
    for entry in range(OG_spectra.shape[1]):
         Match = np.where(np.amin(errorTable) == errorTable)
         matchTable += [Match]
         #print(Match, errorTable[Match])
         errorTable[Match[0],:]=10**7
    
            
    return(matchTable)
def IrPlotter(item,title):
    plt.plot(np.linspace(0,1000,len(item)),item,markersize=.1)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf()
    
    
    
fileList = glob.glob('Tri_A1*/Tri_A1*/input.log')

print(isinstance(fileList[0],str))
plt.vlines(file2Spectra(fileList[1])[:,0],0,file2Spectra(fileList[0])[:,1])


def nmf2TesterMix():
   # pdb.set_trace()

    fraction1 = random.random()
    fraction2 = random.random()
    IR0 =random.choice(fileList)
    IR0 = file2Spectra(IR0)
    print(IR0.shape)
    
    numPoints = len(IR0[:,0])
    IR1 =random.choice(fileList)
    print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}')
    IRF = np.zeros((2,numPoints))
    IROG = np.array([IR0,IR1],dtype='object')
    print("IROG shape", IROG.shape)
    #IrPlotter(IR1,"SecondPlot")
    IRF[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRF[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    
    
    IrMix = np.zeros((4,numPoints))
    IrMix[0,:]=IRF[0,:]
    IrMix[1,:]=IRF[1,:]
   
    #IrMix[3,:] = IR0*fraction1 + IR1*(1-fraction1)
    
    IRF= np.transpose(IRF)
    #IrMix[3,:] = IR0*(1-fraction2)  + IR1*fraction1
   # IrPlotter( IrMix[0,:],"FirstMix")
   # IrPlotter(IrMix[1,:],"SecondMix")
    IrMix= np.transpose(IrMix)
    #model  = NMF(n_components=2, init='nndsvda', max_iter=1000, tol= 1*10**-6, solver='mu')
    #print(model)
    #
   # Wbaby = model.fit_transform(IrMix)
    #Hbaby = model.components_
    model = NMF(n_components=2, max_iter=10000, tol= 1*10**-8, solver= 'cd')#, init='custom')
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    H = model.components_
    HO = H.copy()
    print ("H", HO)
    
    H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    print ("H adjusted", H)
    #print(H)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
    #print(model.fit(IrMix))
    W2 = np.matmul(W,H)

# =============================================================================
    for entri in nmfMatcher(IROG,W):
        
         print('entri', entri)
         plt.plot(np.linspace(0,1000,numPoints),IROG[entri[0][0],:],color="red")
         if H[0,0]>.01:
             print(f'The calculated fraction of the first is {H[0,2]:.5}.')
             print(f'The calculated fraction of the second is {H[1,2]:.5}.')
             

             
         else:
             print(f'The calculated fraction of the first is {H[1,2]:.5}')
             print(f'The calculated fraction of the second is, {H[0,2]:.5}.')
             
         plt.plot(np.linspace(0,1000,numPoints),(W[:,entri[1][0]])*(max(HO[entri[1][0]])))
        # print("full", (max(HO[entri[0][0]])))
        


         plt.gca().invert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()
         
nmf2TesterMix()
