import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

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

def IrPlotter(item,title,ran1,ran2,leg = [], multiple = False):
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,len(item)),item,markersize=.1)
    else:
        for n in item:
            plt.plot(np.linspace(ran1,ran2,len(n)),n,markersize=.1)
    if len (leg) > 0:
        plt.legend(leg,fontsize='x-small')
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")    
    plt.show()
    plt.clf() 
    
    
#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
def gaussian_broadening(spectra, broaden, ran1,ran2,resolution=1):
 
    """ Performs gaussian broadening on IR spectrum
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """
    IR = np.zeros(resolution*(int(ran2-ran1) + 1))
    X = np.linspace(ran1,ran2, resolution*(int(ran2-ran1)+1))

    #for f, i in zip(spectra[:,0]):
      #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
      #  IR=np.vstack((X, IR)).T #tspec
   
    freq = spectra[0]
    inten = spectra[1]
    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
    return IR



#Untreated =gaussian_broadening(fetchIr('UntreatedSample.txt',8),25,1000,4000)
#IrPlotter(Untreated,"untreated",1000,4000)

#print(isinstance(fileList[0],str))
#plt.vlines(file2Spectra(fileList[1])[:,0],0,file2Spectra(fileList[0])[:,1])

def peakPlotter(W,numPeaks,ran1,ran2,col):
    plt.plot(np.linspace(ran1,ran2,(ran2-ran1+1)), W[:,col])
    plt.vlines(numPeaks, 0,W[:,0][numPeaks-ran1]+2)
    plt.title("peaks")
    plt.show()
    plt.clf()



def nmf2TesterMixB(ran1,ran2,broad,res):
    #pdb.set_trace()

    fraction1 = 1
    fraction2 = 1
    Humdities = [5,10,20,20,30,40,50,60,70,80,90,95]
    FileSelection1 = random.randrange(1,12)
    print(f'The first selection is {Humdities[FileSelection1]}% humidity')
    FileSelection2 = 0
    while (FileSelection2 == FileSelection1) or (FileSelection2 == 0) :
        FileSelection2 = random.randrange(1,12)
   
    print(f'The second selection is {Humdities[FileSelection2]}% humidity.')
    #print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}.')
    
    #Creating Two Spectra
    Untreated =gaussian_broadening(fetchIr('UntreatedSample.txt',FileSelection1),broad,ran1,ran2,res)
    IrPlotter(Untreated, f'Untreated Spectra:  {Humdities[FileSelection1]}% humidity',ran1,ran2)

    MEOH1 =gaussian_broadening(fetchIr('MeOHSample.txt',FileSelection1),broad,ran1,ran2,res)
    IrPlotter(MEOH1,f'MeOH Spectra 1:  {Humdities[FileSelection1]}% humidity',ran1,ran2)
    
    MEOH2 =gaussian_broadening(fetchIr('MeOHSample.txt',FileSelection2),broad,ran1,ran2,res)
    IrPlotter(MEOH2,f'MeOH Spectra 2: {Humdities[FileSelection2]}% humidity',ran1,ran2)
    
    WA45_1 =gaussian_broadening(fetchIr('WA45Sample.txt',FileSelection1),broad,ran1,ran2,res)
    IrPlotter(WA45_1,f'WA45 Spectra 1: {Humdities[FileSelection1]}% humidity',ran1,ran2)
    
    WA45_2 =gaussian_broadening(fetchIr('WA45Sample.txt',FileSelection2),broad,ran1,ran2,res)
    IrPlotter(WA45_2,f'WA45 Spectra 2: {Humdities[FileSelection2]}% humidity',ran1,ran2)
    
    
    IrPlotter([Untreated,MEOH1,MEOH2,WA45_1,WA45_2], 'Input Spectra', ran1,ran2,[f'Untreated Spectra:  {Humdities[FileSelection1]}% humidity', 
                                                                                 f'MeOH 1:  {Humdities[FileSelection1]}% humidity',
                                                                                 f'MeOH 2: {Humdities[FileSelection2]}% humidity'
                                                                                 f'WA45 1: {Humdities[FileSelection1]}% humidity'
                                                                                 f'WA45 2: {Humdities[FileSelection2]}% humidity'
                                                                                 ],True)
   # print(IR1.shape)
 
    IRF = np.zeros((5,(ran2-ran1+1)*res))
 #   IRF[0,:] = IR0 *fraction1 + IR1*(1-fraction1)
  #  IRF[1,:] = IR0 * fraction2 +  IR1*(1-fraction2)
   # IrPlotter(IRF[0,:],"fractional Spectra 1")
    #IrPlotter(IRF[1,:],"fractional Spectra 2")
    IRF[0,:] = Untreated
    IRF[1,:] = MEOH1
    IRF[2,:] = MEOH2
    IRF[3,:] = WA45_1
    IRF[4,:] = WA45_2
    
    
    IRF= np.transpose(IRF)
    
    model = NMF(n_components=4, max_iter=5000, tol= 1*10**-10, solver= 'mu', init='nndsvda', beta_loss= 'kullback-leibler')#, alpha = .3  )
    W = model.fit_transform(IRF)
    #IrPlotter(W[:,0])
    
    print ("W-size",W.shape)
    #print("mean", np.mean(W), np.mean(IRF[:,0]))
    numPeaks0 = (find_peaks(W[:,0])[0])+ran1
    
    numPeaks1 = (find_peaks(W[:,1])[0])+ran1
    numPeaks2 = (find_peaks(W[:,2])[0]) +ran1
    numPeaks3 = (find_peaks(W[:,3])[0]) +ran1
    
    peakPlotter(W,numPeaks0,ran1,ran2,0)
    peakPlotter(W,numPeaks1,ran1,ran2,1)
    peakPlotter(W,numPeaks2,ran1,ran2,2)
    peakPlotter(W,numPeaks3,ran1,ran2,3)
    
    
    
    print ("Peaks",numPeaks0,numPeaks1,numPeaks2,numPeaks3 )
    print('nums',numPeaks0,numPeaks1,numPeaks2,numPeaks3)
    K0= np.mean([W]) /  np.mean(IRF[:,0])
    K1= np.mean([W]) /  np.mean(IRF[:,1])
    K2= np.mean([W]) /  np.mean(IRF[:,2])
    K3= np.mean([W]) /  np.mean(IRF[:,3])
    K4= np.mean([W]) /  np.mean(IRF[:,4])
   
     
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,0]*K0 ],'Output Spectra vs Untreated Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "Untreated Spectra"], True)
    
    IrPlotter(W[:,0], "Calculated Spectra 1",ran1,ran2)
    IrPlotter(W[:,1], "Calculated Spectra 2",ran1,ran2)
    IrPlotter(W[:,2], "Calculated Spectra 3",ran1,ran2)
    IrPlotter(W[:,3], "Calculated Spectra 4",ran1,ran2)
    #IrPlotter(W[:,4], "Calculated Spectra 5",ran1,ran2)
    
    
    
    
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,1]*K1 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "MeOH A"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,2]*K2 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "MeOH B"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,3]*K3 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "WA45 A"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,4]*K4 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "WA45 B"], True)
    #IrPlotter([product[:,0],product[:,1],product[:,2],product[:,3],IRF[:,1] ],'4Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", f'MeOH Spectra 1:  {Humdities[FileSelection1]}% humidity'], True)
   # IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3] ],'4Spectra',["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra"], True)
    H = model.components_
    Hnorm = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    print(H)
    print(Hnorm)
    #product = np.matmul(W,H)
    #IrPlotter([product[:,0],product[:,1],product[:,2],product[:,3],IRF[:,0] ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "Untreated Spectra"], True)
    #IrPlotter([product[:,0],product[:,1],product[:,2],product[:,3],IRF[:,1] ],'4Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", f'MeOH Spectra 1:  {Humdities[FileSelection1]}% humidity'], True)
   # matchTable = nmfMatcher (IRF, product)
    #IrOrgs = [IR0,IR1]
    #IrPlotter ([IR0, W[:,0]], "First Calculated Spectra", ['Original Ir', 'NMF Generated IR'],True)
    #IrPlotter ([IR0, W[:,1]], "Second Calculated Spectra", ['Original Ir', 'NMF Generated IR'],True)
   # IrPlotter ([IR0, W[:,2]], "Third Calculated Spectra", ['Original Ir', 'NMF Generated IR'],True)
    #IrPlotter ([IR0, W[:,3]], "Fourth Calculated Spectra", ['Original Ir', 'NMF Generated IR'],True)
    
    #IrPlotter([IrOrgs[matchTable[0,0]], product[:,matchTable[0,1]]],"First Calculated Spectra", ['Original Ir', 'NMF Generated IR'],True)
    #IrPlotter([IrOrgs[matchTable[1,0]] ,product[:,matchTable[1,1]]],"Second Calculated Spectra", ['Original Ir', 'NMF Generated IR'], True)      
   # IrPlotter([IROrgs[matchTable[0][0]],product
               #[W[matchTable[0][1],:]]], "First Matched Spectra", True)
    
    #print(matchTable)
   # 
    
nmf2TesterMixB(1600,1700,20,1)  
   