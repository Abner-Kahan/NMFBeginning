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

def nmfMatcher(OG_spectra,Calc_spectra):
    #print(OG_spectra[:50,0])
    #print("--------------")
    #print(Calc_spectra[650:700,1])
    errorTable = np.zeros((OG_spectra.shape[1], Calc_spectra.shape[1]))
    for n in range (len(OG_spectra)):
        for p in range(OG_spectra.shape[1]):
            for q in range(OG_spectra.shape[1]):
                errorTable[q,p] += abs( OG_spectra[n,q] - Calc_spectra[n,p])
    matchTable=[]
    #print("errorTable \n \n",errorTable)
    for entry in range(OG_spectra.shape[1]):
        Match = np.where(np.amin(errorTable) == errorTable)
        matchTable += [Match]
        #print(Match, errorTable[Match])
        errorTable[Match[0],:]=10**5
    return(matchTable)

def wholeNum(L,multiplier):
#multiplier should be a lower bound int guess
    delta=0
    for n in L:  
        delta  += (n*multiplier - np.floor(n*multiplier))
    if delta > .000001:
        return wholeNum(L,multiplier+1)
    else:
        return multiplier
    
#print(wholeNum([.20,.90,.15],1)  )  
def nmfTester(spectra):
    IR_stack = addIRS(10,20000)
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title("0 of " + str(spectra)+ " : Original Spectra")
    #plt.savefig("0of" + str(spectra)+ ":import matplotlib.pyplot as plt

    plt.show()
    plt.close()
    for n in range(spectra-1):
        IR_stack = np.column_stack((IR_stack,addIRS(10,20000))) 
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,n+1])
        plt.gca().invert_yaxis()
        plt.title(str(n)+ " of " + str(spectra)+ ": Original Spectra")
        #plt.savefig((str(n+1)+ "of" + str(spectra)+ "_ Original Spectra.png"))
        plt.show()
        plt.close()
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Original Spectra")
    plt.xlabel("cm^-1")
   # plt.savefig((str(spectra)+ " Original Spectra.png"))
    plt.show()
    model = NMF(n_components=spectra, init='random', random_state=0, shuffle=1 )
    
    W = model.fit_transform(IR_stack)
 
    for n in range(spectra):
        plt.plot(np.linspace(0,1000,20000),W[:,n],markersize=1)
        plt.gca().invert_yaxis()
        plt.title(str(n)+ " of " + str(spectra)+ ": Calculated Spectra")
       # plt.savefig((str(n)+ " of " + str(spectra)+ "_ Calculated Spectra.png"))
        plt.show()        
        plt.close()
    plt.plot(np.linspace(0,1000,20000),W)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Calculated Spectra")
    plt.xlabel("cm^-1")
    plt.close()
    for entri in nmfMatcher(IR_stack,W):
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,entri[0][0]],color="red")
        plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
        plt.legend(["Original", "Calculated"])
        plt.gca().invert_yaxis()
        plt.title(str(entri[0][0])+ " Both Spectra")
        plt.xlabel("cm^-1")
        plt.show()
        plt.close()
    print(IR_stack.shape)
#nmfTester(3)        
        
def nmfListTester(spectra):
    #IR_stack = addIRS(10,20000)
    #plt.plot(np.linspace(0,1000,20000),IR_stack)
    #plt.gca().invert_yaxis()
    #plt.title("0 of " + str(spectra)+ " : Original Spectra")
    #plt.savefig("0of" + str(spectra)+ ":Original Spectra.png")
   # plt.show()
   
    ind = 0
    for n in spectra:
        n = int(n * wholeNum(spectra, 1))
        spectra[ind] = n
        ind +=1
        
    IR_stack = np.expand_dims(addIRS(1,20),axis=0)
    print (IR_stack)
    print("stop1")

    
    ind= 0
    start = 1
    print (spectra, "spectra")
    for n in spectra:
        copi = np.expand_dims(addIRS(1,20),axis=0)
        print (copi)
        print("stop2")
        if start == 1:  
            for q in range(n-1):
                IR_stack = np.append(IR_stack,copi,axis=0)
            start = 0
            print (IR_stack)
            print("stop3")
        else:
             for q in range(n):
                IR_stack = np.append(IR_stack,copi,axis=0)     
        ind +=1
        print(IR_stack)
        print (IR_stack.shape)
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,n-1])
        plt.gca().invert_yaxis()
        plt.title(str(ind)+ " of " + " Original Spectra")
        #plt.savefig((str(n+1)+ "of" + str(spectra)+ "_ Original Spectra.png"))
        plt.show()
        plt.close()

    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Original Spectra")
    plt.xlabel("cm^-1")
   # plt.savefig((str(spectra)+ " Original Spectra.png"))
  #  plt.show()
  #  model = NMF(n_components=len(spectra), init='random', random_state=0, shuffle=1 )
    '''
    W = model.fit_transform(IR_stack)
    for n in range(len(spectra)):
        plt.plot(np.linspace(0,1000,20000),W[:,n],markersize=1)
        plt.gca().invert_yaxis()
        plt.title(str(n)+ " of "+ ": Calculated Spectra")
       # plt.savefig((str(n)+ " of " + str(spectra)+ "_ Calculated Spectra.png"))
        plt.show()        
        plt.close()
    plt.plot(np.linspace(0,1000,20000),W)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Calculated Spectra")
    plt.xlabel("cm^-1")
    plt.close()
    for entri in nmfMatcher(IR_stack,W):
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,entri[0][0]],color="red")
        plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
        plt.legend(["Original", "Calculated"])
        plt.gca().invert_yaxis()
        plt.title(str(entri[0][0])+ " Both Spectra")
        plt.xlabel("cm^-1")
        plt.show()
        plt.close()
'''
#nmfTester(3)
nmfListTester([.25,.5])
    

    
    #print ("zero-dif",zz)
   # print("one-dif",zo)
    #print(OG_spectra.shape[1])
    #plt.plot(np.linspace(0,1000,20000),OG_spectra)
    #plt.plot(np.linspace(0,1000,20000),Calc_spectra)
#nmfMatcher(outboy[0],outboy[1])    
    
#Changes
#Number of peaks changed to 10
# Changed to Gaussian
