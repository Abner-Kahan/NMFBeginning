#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF
#import seaborn as sns
#import matplotlib.ticker as mtick
import scipy.integrate as integrate
import scipy.special as special


def addIRS(peakNum,PointNum):
    Ir = np.zeros(PointNum)
    Irlocs = np.zeros(PointNum)

    for n in range(peakNum):
        thickness = int (PointNum * random.choice([.007,.008,.009,.01,.011,.012,.013,.014,.042]))
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
       # Ir[start:end] *= 1/ np.max(Ir)
        plt.clf()
    return Ir
#IR_stack = addIRS(10,10000)
#plt.gca().invert_yaxis()
#plt.plot(np.linspace(0,1000,10000),IR_stack)

#IR2Stack= addIRS(30,20000).append[(addIRS(30,20000))]
#IR2Stack = np.column_stack((addIRS(30,20000),addIRS(30,20000)))

def nmfMatcher(OG_spectra,Calc_spectra):
    #Sprint(OG_spectra[:10,:])
    True
    #print(Calc_spectra[:10,:])
    
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
    #print("bob \n", matchTable)
    #for match in matchTable:
      #if matchTable.count(match) > 1:
        #  matchTable.remove(match)
    
            
    return(matchTable)

#print(wholeNum([.20,.90,.15],1)  )  
def nmfTester(spectra):
    IR_stack = addIRS(10,20000)
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title("0 of " + str(spectra)+ #sns.distplot(IR_stack,bins=1000)
" : Original Spectra")
    #plt.savefig("0of" + str(spectra)+ ":import matplotlib.pyplot as plt

    plt.show()
    plt.clf()
    for n in range(spectra-1):
        IR_stack = np.column_stack((IR_stack,addIRS(10,20000))) 
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,n+1])
        plt.gca().invert_yaxis()
        plt.title(str(n+1)+ " of " + str(spectra)+ ": Original Spectra")
        #plt.savefig((str(n+1)+ "of" + str(spectra)+ "_ Original Spectra.png"))
        plt.show()
        plt.clf()
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Original Spectra")
    plt.xlabel("cm^-1")
   # plt.savefig((str(spectra)+ " Original Spectra.png"))
    plt.show()
    model = NMF(n_components=spectra, init='random', random_state=0, shuffle=1 )
    #print(IR_stack[:50,:])
    W = model.fit_transform(IR_stack)
 
    for n in range(spectra):
        plt.plot(np.linspace(0,1000,20000),W[:,n],markersize=1)
        plt.gca().invert_yaxis()
        plt.title(str(n)+ " of " + str(spectra)+ ": Calculated Spectra")
       # plt.savefig((str(n)+ " of " + str(spectra)+ "_ Calculated Spectra.png"))
        plt.show()        
        plt.clf()
    plt.plot(np.linspace(0,1000,20000),W)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Calculated Spectra")
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf()
    for entri in nmfMatcher(IR_stack,W):
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,entri[0][0]],color="red")
        plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
        plt.legend(["Original", "Calculated"])
        plt.gca().invert_yaxis()
        plt.title(str(entri[0][0])+ " Both Spectra")
        plt.xlabel("cm^-1")
        plt.show()
        plt.clf()
    print(IR_stack.shape)
#nmfTester(3)        
def IrPlotter(item,title):
    plt.plot(np.linspace(0,1000,len(item)),item,markersize=.2)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf()
    
def nmfTester2Frac(fraction):
    IrMix = np.zeros((2,20000))
    IrMix[0,:] = addIRS(10,20000)*fraction
    IrMix[1,:] = addIRS(10,20000)*(1-fraction)
    IrPlotter(IrMix[0,:], "IrMixSpectra Small")
    IrPlotter(IrMix[1,:], "IrMixSpectra Large")
    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd', random_state=0, max_iter=6000 )
    W = model.fit_transform(IrMix)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0))))
    #print(model.fit(IrMix))
    for entri in nmfMatcher(IrMix,W):
        
        plt.plot(np.linspace(0,1000,20000),IrMix[:,entri[0][0]],color="red")
        plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
        plt.gca().invert_yaxis()
        plt.legend(["Original", "Calculated"])
        plt.title(str(entri[0][0])+ " Both Spectra")
        plt.xlabel("cm^-1")
        plt.show()
        plt.clf()
        #print( np.average(W[:,entri[1][0]]) / np.average( W[:,entri[0][0]]))
     
    
#nmfTester2Frac(.10)

    

def nmfTesterFrac(fractionList):
     leng = len(fractionList)
    
     IrMix = np.zeros((leng,20000))
     
     for n in range(leng):
         IrMix[n,:] = addIRS(10,20000)#* fractionList[n]
         IrPlotter(IrMix[n,:], "IrMixSpectra" + str(n) )
     IrMix= np.transpose(IrMix)
     model = NMF(n_components=leng, init='nndsvd', random_state=0, max_iter=500000 )
     W = model.fit_transform(IrMix)
     for n in range(leng):
         IrPlotter(W[:,n], "Calc Spectra"+  str(n))

     for entri in nmfMatcher(IrMix,W):
        
         plt.plot(np.linspace(0,1000,20000),IrMixTrue[:,entri[0][0]],color="red")
         plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
         plt.gca().invert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()
         
         
        
# =============================================================================
#=============================================================================
#nmfTesterFrac([.3,.4,.8,.54,.32,.10,.68,.10,.03,.29,.98,.02,.92,.04])

def nmf2TesterMix(fraction1,fraction2,randi):
    IRF = np.zeros((2,20000))
    IR0 = addIRS(10,20000)
    IrPlotter(IR0,"FirstPlot")
    IR1 = addIRS(10,20000)
    IrPlotter(IR1,"SecondPlot")
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((2,20000))
    IrMix[0,:] = IR0*fraction1 + IR1*(1-fraction1) 
    IrMix[1,:] = IR0*(1-fraction2)  + IR1*fraction1
    IrPlotter( IrMix[0,:],"FirstMix")
    IrPlotter(IrMix[1,:],"SecondMix")
    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd',  max_iter=4000, tol= 1*10**-6, solver='cd' )
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
    #print(model.fit(IrMix))
# =============================================================================
    for entri in nmfMatcher(IRF,W):
        
         plt.plot(np.linspace(0,1000,20000),IRF[:,entri[0][0]],color="red")
         plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
         plt.gca().invert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()
# =============================================================================
nmf2TesterMix(.2,.82,False)

def nmfSimplerMix(fraction):
    IRF = np.zeros((2,20000))
    IR0 = addIRS(10,20000)
    IrPlotter(IR0,"FirstPlot")
    IR1 = addIRS(10,20000)
    IrPlotter(IR1,"SecondPlot")
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((2,20000))
    IrMix[0,:] = IR0*fraction + IR1*(1-fraction) 
    IrMix[1,:] = (IrMix[0,:])+.01
    IrPlotter( IrMix[0,:],"Mix")
    
    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd', max_iter=1000 )
    W = model.fit_transform(IrMix)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0))))
    #print(model.fit(IrMix))
# =============================================================================
    for entri in nmfMatcher(IRF,W):
        
         plt.plot(np.linspace(0,1000,20000),IRF[:,entri[0][0]],color="red")
         plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
         plt.gca().invert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()
#nmfSimplerMix(.20)
#IR comp1
#IR comp2
#mix 1
#mix 2
#GLYCLOPIDIS - FIGURE
#nmfListTester([.25,.75])
    


    
#Changes
#Number of peaks changed to 10
# Changed to Gaussian
