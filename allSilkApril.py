# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 08:20:12 2022

@author: Abner
"""


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


amide1_ran1 =1590
amide1_ran2 =1710


Humdities = [5,10,20,30,40,50,60,70,80,90,95]
solvents = ["Untreated", "MeOH", "WA45"]
broad = 15
x_range = np.linspace(amide1_ran1,amide1_ran2,(amide1_ran2-amide1_ran1+1))

def fetchIr(path,column,ran1,ran2):
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
            if (x == 0) and (float(word) < ran1 or float(word) > ran2):
                
                break


            #print(word, '\n')
            #if float(word) < ran1:
          #      break
           # if float(word) > ran2:
            #    break
            IrArray[x,y] = word
            x+=1
        y+=1
        x=0
    #print(IrArray.shape)
    mask = (IrArray[0,:]) != [0] *len(IrArray[0])
    #print(mask)
    
    IrArray2 = IrArray[:,mask]
    #print(IrArray2.shape)
    #plt.scatter(IrArray2[0],IrArray2[1])
    #print(IrArray)
    
    
    return IrArray2[(0,column),:]

#print(fetchIr('UntreatedSample.txt',3))

def IrPlotter(item,title,ran1,ran2, leg = [], multiple = False):
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
def gaussian_broadening(spectra, broaden, ran1,ran2,experiment =False, resolution=1):
 
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
    if experiment:
         freq = spectra[0] *.965
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
    plt.vlines(numPeaks, 0,W[:,col][numPeaks-ran1])
    plt.legend(np.ndarray.tolist(numPeaks))
    plt.ylim(top =2)
    plt.title(f"peaks Calculated Spectra{col}")
    plt.show()
    plt.clf()


sumsol  = np.load('sumsol.npy')


                  
                  

def nmf2TesterMixB(broad):
    #pdb.set_trace()
    IRF = np.zeros((33,(amide1_ran2-amide1_ran1+1)))
    #change this
    #IRF = np.zeros((22,(amide1_ran2-amide1_ran1+1)))
    #IRF = np.zeros((11,(amide1_ran2-amide1_ran1+1)))
    for n in range(11):
            #change these
           IRF [n,:]= gaussian_broadening(fetchIr('UntreatedSample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)
           #IRF [n,:]= gaussian_broadening(fetchIr('MeOHSample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)
          # IRF [n,:]= gaussian_broadening(fetchIr('WA45Sample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)
           
           #IRF[n+11,:] =  gaussian_broadening(fetchIr('MeOHSample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)
           IRF[n+11,:] =  gaussian_broadening(fetchIr('WA45Sample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)
           
           IRF[n+22,:] =  gaussian_broadening(fetchIr('WA45Sample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)

    Humdities = [5,10,20,30,40,50,60,70,80,90,95]
    fig, axs = plt.subplots(4)
    
    figB, axsB = plt.subplots(nrows=1, ncols=2)
 
    for spec in IRF[:11,:]:
        axs[0].plot(x_range,spec)
        
    axs[0].set_title("Untreated Sample", x=.5, y= .2)
    
    #axs[0].legend(Humdities)

    
    for spec in IRF[11:22,:]:   
        axs[1].plot(x_range, spec)
    
    axs[1].set_title("MeOH Sample", x=.5, y= .2)
    
    #axs[1].legend(Humdities)

    for spec in IRF[22:33,:]:   
        axs[2].plot(x_range, spec)
    
    axs[2].set_title("WA45ample Sample", x=.5, y= .2)
    
    #axs[2].legend(Humdities)
    randSpec = random.choice(list(range(11))+list( range(22,33)))
    axsB[0].plot(IRF[randSpec,:])

    
    #IrPlotter( IRF [:11,:], 'Unstreated Spectr IRF[n+22,:] =  gaussian_broadening(fetchIr('WA45Sample.txt',n+1,amide1_ran1,amide1_ran2),broad,amide1_ran1,amide1_ran2)a', amide1_ran1,amide1_ran2, ['5','10','20','30','40','50','60','70','80','90','95']
           #                                                                     ,True)
    #IrPlotter( IRF [11:22,:], 'MEOH Spectra', amide1_ran1,amide1_ran2, ['5','10','20','30','40','50','60','70','80','90','95']
       #                                                                          ,True)
    #IrPlotter( IRF [22:,:], 'WA45 Spectra', amide1_ran1,amide1_ran2,['5','10','20','30','40','50','60','70','80','90','95']
       #                                                                          ,True)
   # print(IR1.shape)
    

   # IrPlotter (gaussian_broadening(fetchIr('UntreatedSample.txt',1),broad,ran1,ran2,res), "test", ran1, ran2)
    IRF= np.transpose(IRF)
    #model = NMF(n_components=4, max_iter=3000, tol= 1*10**-12, solver= 'cd', init= "nndsvdar", beta_loss= 'frobenius')
    model = NMF(n_components=4, max_iter=3000, tol= 1*10**-12, solver= 'mu', init='nndsvda', beta_loss= 'kullback-leibler')#, alpha = .3  )
    W = model.fit_transform(IRF)
    #IrPlotter(W[:,0])
    
    print ("W-size",W.shape)
    #print("mean", np.mean(W), np.mean(IRF[:,0]))
    numPeaks0 = (find_peaks(W[:,0],'prominence' == 1)[0])+amide1_ran1
    
    numPeaks1 = (find_peaks(W[:,1],'prominence' == 1)[0])+amide1_ran1
    numPeaks2 = (find_peaks(W[:,2],'prominence' == 1)[0]) +amide1_ran1
    numPeaks3 = (find_peaks(W[:,3],'prominence' == 1)[0]) +amide1_ran1
    #numPeaks4 = (find_peaks(W[:,4],'prominence' == .1)[0]) +amide1_ran1
    
# =============================================================================
#     peakPlotter(W,numPeaks0,amide1_ran1,amide1_ran2,0)
#     peakPlotter(W,numPeaks1,amide1_ran1,amide1_ran2,1)
#     peakPlotter(W,numPeaks2,amide1_ran1,amide1_ran2,2)
#     peakPlotter(W,numPeaks3,amide1_ran1,amide1_ran2,3)
# =============================================================================
    
    axs[3].legend(["NMF Calculated 1", "NMF Calculated 2", "NMF Calculated 3", "NMF Calculated 4"])
    axs[3].set_title("NMF Decomposition",  x=.5, y= .05)
    print ("Peaks",numPeaks0,numPeaks1,numPeaks2,numPeaks3 )
    print('nums',numPeaks0,numPeaks1,numPeaks2,numPeaks3)
    K0= np.mean([W]) /  np.mean(IRF[:,0])
    K1= np.mean([W]) /  np.mean(IRF[:,1])
    K2= np.mean([W]) /  np.mean(IRF[:,2])
    K3= np.mean([W]) /  np.mean(IRF[:,3])
    #K4= np.mean([W]) /  np.mean(IRF[:,4])
   
     
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,0]*K0 ],'Output Spectra vs Untreated Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "Untreated Spectra"], True)
    
    #IrPlotter(W[:,ordPeaks[0]], "NMF Calculated 1",amide1_ran1,amide1_ran2)
    #IrPlotter(W[:,ordPeaks[1]], "NMF Calculated 2",amide1_ran1,amide1_ran2)
    #IrPlotter(W[:,ordPeaks[2]], "NMF Calculated 3",amide1_ran1,amide1_ran2)
    #IrPlotter(W[:,ordPeaks[3]], "NMF Calculated 4",amide1_ran1,amide1_ran2)
    #IrPlotter(W[:,4], "NMF Calculated 5",amide1_ran1,amide1_ran2)
    
    #IrPlotter(W[:,4], "Calculated Spectra 5",ran1,ran2)
    SingPeaks =[]
    SingPeaks += [np.argmax(W[:,0])]
    SingPeaks += [np.argmax(W[:,1])]
    SingPeaks += [np.argmax(W[:,2])]
    SingPeaks += [np.argmax(W[:,3])]
    print(SingPeaks, "\n\n\n\n\n----\n\n---")
    sorter =np.argsort(SingPeaks)
    
    W_test = W
    print('\n\n\n\n' ,W, '')
    axs[3].plot(x_range, W_test[:,sorter[0]])
    #meany = np.mean(W_test[:,sorter[0]]/sumsol[0])

    #plt.plot(x_range, sumsol[0]*meany)

    
    axs[3].plot(x_range, W_test[:,sorter[1]])
    #meany = np.mean(W_test[:,sorter[1]]/sumsol[1])

    #plt.plot(x_range, sumsol[1]*meany)

    
    axs[3].plot(x_range, W_test[:,sorter[2]])
   # meany = np.mean(W_test[:,sorter[2]]/sumsol[2])

   # plt.plot(x_range, sumsol[2]*meany)

    
    axs[3].plot(x_range, W_test[:,sorter[3]])
   # meany = np.mean(W_test[:,sorter[3]]/sumsol[3])
    axs[0].set_xticklabels(())
    axs[1].set_xticklabels(())
    axs[2].set_xticklabels(())
    #axs[0].set(xticks=None)
    #axs[1].set(xticks=None)
    #axs[2].set(xticks=None)
   # plt.plot(x_range, sumsol[3]*meany)
    axs[3].set(xlabel='cm$^-1$')
    fig.suptitle("NMF Decomposition of Amide II of Silk")
    #plt.legend(["NMF", "Gaussian"])

    #IrPlotter(sumsol[0], "Gaussian",amide1_ran1,amide1_ran2) 
    axsB[0].plot(W_test[:,sorter[0]])
    axsB[0].plot(W_test[:,sorter[1]])
    axsB[0].plot(W_test[:,sorter[2]])
    axsB[0].plot(W_test[:,sorter[3]])
    
    
    
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,1]*K1 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "MeOH A"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,2]*K2 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "MeOH B"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,3]*K3 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "WA45 A"], True)
    #IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3],IRF[:,4]*K4 ],'Output Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", "WA45 B"], True)
    #IrPlotter([product[:,0],product[:,1],product[:,2],product[:,3],IRF[:,1] ],'4Spectra',ran1,ran2,["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra", f'MeOH Spectra 1:  {Humdities[FileSelection1]}% humidity'], True)
   # IrPlotter([W[:,0],W[:,1],W[:,2],W[:,3] ],'4Spectra',["1st Spectra", "2nd Spectra", "3rd Spectra", "4th Spectra"], True)
    H = model.components_
    axsB[1].plot(np.matmul(W,H))
    plt.show()
    #plt.plot(np.dot(W,H[0,randSpec]))
    #Hnorm = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    # for n in range(11): 
    #     plt.plot(np.dot(W,H[0,n]))
    #     #plt.plot(np.dot(W,H[1,n]))
    #     #plt.plot(np.dot(W,H[2,n]))
    #    # plt.plot(np.dot(W,H[3,n]))
        
    # plt.title('NMF Rebuilt Untreated Sample')
    # plt.xlabel("cm-1")
    # plt.legend(Humdities)
    # plt.show()
    # plt.clf()
    
    # for n in range(11,22): 
    #     plt.plot(np.dot(W,H[0,n]))
    # plt.title('NMF Rebuilt MeOH Sample')
    # plt.xlabel("cm-1")
    # plt.legend(Humdities)
    # plt.show()
    # plt.clf()
    
    # for n in range(22,33): 
    #     plt.plot(np.dot(W,H[0,n]))
    # plt.title('NMF Rebuilt WA45 Sample')
    # plt.xlabel("cm-1")
    # plt.legend(Humdities)
    # plt.show()
    # plt.clf()
    #np.save("H.npy",H)
    #np.save("Hnorm.npy",Hnorm)
   #  Hnorm = np.apply_along_axis(lambda l : l/np.sum(l),0,H)
   #  plt.plot(Humdities, Hnorm[sorter[0], :11], label = numPeaks0)
   #  plt.plot(Humdities, Hnorm[sorter[1], :11], label = numPeaks1)
   #  plt.plot(Humdities, Hnorm[sorter[2], :11], label = numPeaks2)
   #  plt.plot(Humdities, Hnorm[sorter[3], :11], label = numPeaks3)
   # # plt.plot(Humdities, Hnorm[4, :11], label = numPeaks4)
   #  plt.xlabel("%humidity")
   #  plt.legend()
   #  plt.title('Untreated Sample')
   #  plt.show()
   #  plt.clf()
    
   #  plt.plot(Humdities, Hnorm[sorter[0], 11:22], label = numPeaks0)
   #  plt.plot(Humdities, Hnorm[sorter[1], 11:22], label = numPeaks1)
   #  plt.plot(Humdities, Hnorm[sorter[2], 11:22], label = numPeaks2)
   #  plt.plot(Humdities, Hnorm[sorter[3], 11:22], label = numPeaks3)
   # # plt.plot(Humdities, Hnorm[4, 11:22], label = numPeaks4)

   #  plt.xlabel("%humidity")
   #  plt.legend()
   #  plt.title('MeOH Fractions')
   #  plt.show()
   #  plt.clf()    
    
   #  plt.plot(Humdities, Hnorm[sorter[0], 22:33], label = numPeaks0)
   #  plt.plot(Humdities, Hnorm[sorter[1], 22:33], label = numPeaks1)
   #  plt.plot(Humdities, Hnorm[sorter[2], 22:33], label = numPeaks2)
   #  plt.plot(Humdities, Hnorm[sorter[3], 22:33], label = numPeaks3)
   # # plt.plot(Humdities, Hnorm[4, 22:33], label = numPeaks4)

   #  plt.xlabel("%humidity")
   #  plt.legend()
   #  plt.title('WA45 Fractions')
   #  plt.show()
   #  plt.clf()
    
    
    
    return H

H = nmf2TesterMixB(15)  