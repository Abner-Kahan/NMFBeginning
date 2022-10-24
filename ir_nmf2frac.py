#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF

#np.set_printoptions(precision=3)

def addIRS(peakNum,PointNum):
    #Ir = np.zeros(PointNum)
    #Irlocs = np.zeros(PointNum)

    #for n in range(peakNum):
    #    thickness = int (PointNum * random.choice([.007,.008,.009,.01,.011,.012,.013,.014,.042]))
        #print (thickness)
        #for b in range(thickness):
        #    if any(Ir[b:b+10]):
        #        Irlocs[b] = 1
        #for c in range(PointNum - thickness):
        #    if any(Ir[c + thickness : (c +thickness+10)]):
        #        Irlocs[c] = 1
                
    #    xloc=  random.choice(np.where(Irlocs == 0)[0])
        #print("xloc", xloc)
         #frcations of graph
        #print("thickness", thickness)
        #thickness2 = 1. #random.randrange(1,5)
     #   peakHeight=random.random()
#This is supposed to deal with peaks toward the ends of the spectra
      #  if (xloc + 1 + thickness ) > PointNum:
       #     end = PointNum
       # else: 
       #     end = (xloc + 1 + thickness)
       # if (xloc -  thickness < 0): 
       #     start = 0
       # else:
       #start = xloc - thickness
        #print("start", start)
        #print("end", end)
        #start = int(start)    axs[0,1].set(xticks = [])
        #end = int(end)
        #Ir[start:end]= np.random.normal((end-start)/2,thickness,end-start)
        #Ir[start:end] = stats.norm.pdf(np.linspace(start,end,end-start),(end+start)/2,thickness)*peakHeight
        #mean = int((end+start)/2)
        #print(mean)
        #Ir[mean-thickness:mean+thickness+1] = stats.norm.pdf(np.linspace(start,end,end-start), mean, thickness/3)
# Ir[start:end] *= 1/ np.max(Ir)
        #plt.clf()
        #return Ir


        #IR = np.zeros((int(4000/resolution) + 1,))
        #X = np.linspace(0,4000, int(4000/resolution)+1)
        IR = np.zeros((PointNum,))
        X =  np.linspace(0, PointNum, PointNum)

        xloc = [ random.randrange(0, PointNum) for i in range(peakNum) ] 
        peakHeigh = [ random.random() for i in range(peakNum) ] 
        thickness = [ PointNum * random.choice([.007,.008,.009,.01,.011,.012,.013,.014,.042]) for i in range(peakNum) ] 

        for f, i, b in zip(xloc, peakHeigh, thickness):  
            IR += i*np.exp(-0.5*((X-f)/int(b))**2)

        plt.clf
        return IR 

        #self.IR=np.vstack((X, IR)).T #tspec






def nmfMatcher(OG_spectra,Calc_spectra):
 


    #print(len(OG_spectra))
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
    plt.plot(np.linspace(0,1000,len(item)),item,markersize=.2)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf() 
    
    
def nmf2TesterMix(numPoints):
    fraction1 = random.random()
    fraction2 = random.random()
    print(f'The expected fractions are  {fraction1:.3}, {fraction2:.3}.')
    fig, axs = plt.subplots(3,2)
    IRF = np.zeros((2,numPoints))
    IR0 = addIRS(8,numPoints)
    axs[0,0].plot(IR0/np.max(IR0), color= 'goldenrod')
    axs[0,0].set_title("Glycan 1", x=.5, y= .75)
    #IrPlotter(IR0,"FirstPlot")
    IR1 = addIRS(8,numPoints)
    axs[0,1].plot(IR1/np.max(IR1), color= 'dodgerblue')
    axs[0,1].set_title("Glycan 2", x=.5, y= .75)
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((2,numPoints))
    IrMix[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMix[1,:] = IR0*fraction2 + IR1*(1-fraction2) 
    #IrMix[3,:] = IR0*(1-fraction2)  + IR1*fraction1
    axs[1,0].plot(IrMix[0,:]/np.max(IrMix[0,:]), color= 'springgreen')
    axs[1,0].set_title("Mixture 1", x=.5, y= .75)
    axs[1,1].plot(IrMix[1,:]/np.max(IrMix[1,:]), color= 'springgreen')
    axs[1,1].set_title("Mixture 2", x=.5, y= .75)
    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd',  max_iter=15000, tol= 1*10**-7)
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    H = model.components_
    HO = H.copy()
    print(H)
    H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    #print(H)
    axs[2,0].plot(W[:,0]/np.max(W[:,0]), color= 'goldenrod')
    axs[2,0].set_title("NMF Glycan 1 ", x=.5, y= .75)
    axs[2,1].plot(W[:,1]/np.max(W[:,1]), color= 'dodgerblue')
    axs[2,1].set_title("NMF Glycan 2", x=.5, y= .75)
    axs[0,0].set(xticks = [])
    axs[1,0].set(xticks = [])
    axs[0,1].set(xticks = [])
    axs[1,1].set(xticks = [])
    axs[2,0].set(xlabel='cm$^-1$')
    axs[2,1].set(xlabel='cm$^-1$')
    axs[2,0].set(xticks = [np.linspace(0, numPoints, 11).all()])
    axs[2,1].set(xticks = [np.linspace(0, numPoints, 11).all()])
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
nmf2TesterMix(10000)
