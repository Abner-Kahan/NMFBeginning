import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF

def fetchIr(path,column):
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
            IrArray[x,y] = word
            x+=1
        y+=1
        x=0
    return IrArray[(0,column),:]

#print(fetchIr('UntreatedSample.txt',3))

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
    
    
#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
def gaussian_broadening(spectra, broaden, resolution=1):
 
    """ Performs gaussian broadening on IR spectrum
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """
    IR = np.zeros((int(4000/resolution) + 1))
    X = np.linspace(0,4000, int(4000/resolution)+1)

    #for f, i in zip(spectra[:,0]):
      #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
      #  IR=np.vstack((X, IR)).T #tspec
   
    freq = spectra[0]
    inten = spectra[1]
    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
    return IR


def nmfMatcher(OG_spectra,Calc_spectra):
    pdb.set_trace() 
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




def nmf2TesterMixB():
    #pdb.set_trace()

    fraction1 = random.random()
    fraction2 = random.random()
    FileSelection1 = random.randrange(1,12)
    print(f'The first selection is {FileSelection1}.')
    FileSelection2 = random.randrange(1,12)
    print(f'The second selection is {FileSelection2}.')
    #print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}.')
    
    #Creating Two Spectra
    IR0 =gaussian_broadening(fetchIr('MeOHSample.txt',FileSelection1),25)
    IrPlotter(IR0,"MeOH Spectra")

    IR1 =gaussian_broadening(fetchIr('UntreatedSample.txt',FileSelection2),25)
    IrPlotter(IR1,"Untreated Spectra")
   # print(IR1.shape)
 
    IRF = np.zeros((2,4001))
    IRF[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
    IRF[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
    IrPlotter(IRF[0,:],"fractional Spectra 1")
    IrPlotter(IRF[1,:],"fractional Spectra 2")
    IRF= np.transpose(IRF)
    
    model = NMF(n_components=2, max_iter=4000, tol= 1*10**-10, solver= 'mu', init ='nndsvda', beta_loss= 'kullback-leibler' )
    W = model.fit_transform(IRF)
    H = model.components_
    Hnorm = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    print(H)
    print(Hnorm)
    product = np.matmul(W,H)
    matchTable = nmfMatcher (IRF, product)
    IrOrgs = [IR0,IR1]
    IrPlotter([IrOrgs[matchTable[0,0]], product[:,matchTable[0,1]]],"First Calculated Spectra", ['Original Ir', 'NMF Generated IR'],True)
    IrPlotter([IrOrgs[matchTable[1,0]] ,product[:,matchTable[1,1]]],"Second Calculated Spectra", ['Original Ir', 'NMF Generated IR'], True)      
   # IrPlotter([IROrgs[matchTable[0][0]],product
               #[W[matchTable[0][1],:]]], "First Matched Spectra", True)
    
    print(matchTable)
nmf2TesterMixB()


#IrPlotter(gaussian_broadening(fetchIr('MeOHSample.txt',1),25), "test1")
#IrPlotter(gaussian_broadening(fetchIr('UntreatedSample.txt',2),25), "test1")
#IrPlotter(gaussian_broadening(fetchIr('WA45Sample.txt',3),25), "test1")
#IrPlotter(gaussian_broadening(fetchIr('UntreatedSample.txt',4),25), "test1")
#IrPlotter(gaussian_broadening(fetchIr('UntreatedSample.txt',5),25), "test1")
#IrPlotter(gaussian_broadening(fetchIr('UntreatedSample.txt',6),25), "test1")


