#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF
#import seaborn as sns
import matplotlib.ticker as mtick



def addIRS(peakNum,PointNum):
    Ir = np.zeros(PointNum)
    Irlocs = np.zeros(PointNum)

    for n in range(peakNum):
        thickness = int (PointNum * random.choice([.007,.008,.009,.01,.011,.012,.013,.014,.06]))
        for b in range(thickness):
            if any(Ir[b:b+10]):
                Irlocs[b] = 1
        for c in range(PointNum - thickness):
            if any(Ir[c + thickness : (c +thickness+10)]):
                Irlocs[c] = 1
                

        xloc=  random.choice(np.where(Irlocs == 0)[0])
        #print("xloc", xloc)
         #frcations of graph
        #print("thickness", thickness)
        #thickness2 = 1. #random.randrange(1,5)
        peakHeight=random.random()
#This is supposed to deal with peaks toward the ends of the spectra
        if (xloc + 1 + thickness ) > PointNum:
            end = PointNum
        else: 
            end = (xloc + 1 + thickness)
        if (xloc -  thickness < 0):
            start = 0
        else:
            start = xloc - thickness
        #print("start", start)
        #print("end", end)
        start = int(start)
        end = int(end)
        #Ir[start:end]= np.random.normal((end-start)/2,thickness,end-start)
        Ir[start:end] = stats.norm.pdf(np.linspace(start,end,end-start),(end+start)/2,thickness)*thickness*peakHeight
    return Ir
#IR_stack = addIRS(10,10000)
#plt.gca().invert_yaxis()
#plt.plot(np.linspace(0,1000,10000),IR_stack)

#sns.distplot(IR_stack,bins=1000)
#IR2Stack= addIRS(30,20000).append[(addIRS(30,20000))]
#IR2Stack = np.column_stack((addIRS(30,20000),addIRS(30,20000)))



def nmfTester(spectra):
    IR_stack = addIRS(10,20000)
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title("0 of " + str(spectra)+ " : Original Spectra")
    plt.savefig("0of" + str(spectra)+ ":Original Spectra.png")
    plt.show()
    plt.close()
    for n in range(spectra - 1):
        IR_stack = np.column_stack((IR_stack,addIRS(10,20000)))
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,n+1])
        plt.gca().invert_yaxis()
        plt.title(str(n+1)+ " of " + str(spectra)+ ": Original Spectra")
        plt.savefig((str(n+1)+ " of " + str(spectra)+ "_ Original Spectra.png"))
        plt.show()
        plt.close()
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Original Spectra")
    plt.xlabel("cm^-1")
    plt.savefig((str(spectra)+ " Original Spectra.png"))
    plt.show()
    model = NMF(n_components=spectra, init='random', random_state=0,alpha=1,max_iter=3000 )
    W = model.fit_transform(IR_stack)
    for n in range(spectra):
        plt.plot(np.linspace(0,1000,20000),W[:,n],markersize=1)
        plt.gca().invert_yaxis()
        plt.title(str(n)+ " of " + str(spectra)+ ": Calculated Spectra")
        plt.savefig((str(n)+ " of " + str(spectra)+ "_ Calculated Spectra.png"))
        plt.show()
        plt.close()
    plt.plot(np.linspace(0,1000,20000),W)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Calculated Spectra")
    plt.xlabel("cm^-1")
    return IR_stack,W


outboy = nmfTester(3)

def nmfMatcher(OG_spectra,Calc_spectra):
    #print(OG_spectra[:50,0])
    #print("--------------")
    #print(Calc_spectra[650:700,1])
    errorTable = np.zeros((OG_spectra.shape[1], Calc_spectra.shape[1]))
    for n in range (len(OG_spectra)):
        for p in range(OG_spectra.shape[1]):
            for q in range(OG_spectra.shape[1]):
                errorTable[q,p] += abs(  OG_spectra[n,q] - Calc_spectra[n,p])
    print(errorTable)
    
    #print ("zero-dif",zz)
   # print("one-dif",zo)
    #print(OG_spectra.shape[1])
    #plt.plot(np.linspace(0,1000,20000),OG_spectra)
    #plt.plot(np.linspace(0,1000,20000),Calc_spectra)
       
nmfMatcher(outboy[0],outboy[1] )    
    
    
#Changes
#Number of peaks changed to 10
# Changed to Gaussian
