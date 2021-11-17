#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF

#np.set_printoptions(precision=3)

def addIRS(peakNum,graph):
   

        #IR = np.zeros((int(4000/resolution) + 1,))
        #X = np.linspace(0,4000, int(4000/resolution)+1)
        IR = np.zeros((graph,))
        X =  np.linspace(0, graph, graph)

        xloc = [ random.randrange(0, graph) for i in range(peakNum) ] 
        peakHeight = [ random.random() for i in range(peakNum) ] 
        thickness = [ graph * random.randrange(7,28)/1000 for i in range(peakNum) ] 

        for f, i, b in zip(xloc, peakHeight, thickness):  
            IR += i*np.exp(-0.5*((X-f)/int(b))**2)

        plt.clf
        return IR 

        #self.IR=np.vstack((X, IR)).T #tspec


def nmfMatcher(OG_spectra,Calc_spectra):
 
    
    #print(len(OG_spectra))
    OG_spectra = np.transpose(OG_spectra)
    errorTable = np.zeros((OG_spectra.shape[1], Calc_spectra.shape[1]))
    for n in range (OG_spectra.shape[0]):
         for p in range(OG_spectra.shape[1]):
             for q in range(Calc_spectra.shape[1]):
                 errorTable[p,q] += abs( OG_spectra[n,p] - Calc_spectra[n,q])
   # print("hi \n", errorTable)
    matchTable=[]
    #print("errorTable \n \n",errorTable)
    for entry in range(OG_spectra.shape[1]):
         Match = np.where(np.amin(errorTable) == errorTable)
         matchTable += [Match]
         #print(Match, errorTable[Match])
         errorTable[Match[0],:]=10**7
   # print("Here are matches",matchTable)
            
    return(matchTable)
  
def IrPlotter(item,title):
    plt.plot(np.linspace(0,1000,len(item)),item,markersize=.1)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf()
    
    
def nmf2TesterMix(numPoints,Nfrac):
    

    
    IRMix = np.zeros((Nfrac,numPoints))
    #print(IRMix.shape[0])
   
    
    IR0 = addIRS(8,numPoints)
    IrPlotter(IR0,"First Plot")
    #IrPlotter(IR0,"FirstPlot")
    IR1 = addIRS(8,numPoints)
    IrPlotter(IR1,"Second Plot")

    randtable= []
    for row in range(Nfrac):
         rand1 = random.random()
         rand2 = random.random()
         IRMix[row,:] = IR0*rand1 + IR1 *rand2
         randtable+= [[rand1,rand2]]
         IrPlotter(IRMix[row,:],("Mix"+ str(row)))
    print(randtable)
    IROG = np.array([IR0,IR1])
    #print("IROG shape", IROG.shape)
    
    IRMix= np.transpose(IRMix)
    #model  = NMF(n_components=2, init='nndsvda', max_iter=1000, tol= 1*10**-6, solver='mu')
    #print(model)
    #
   # Wbaby = model.fit_transform(IrMix)
    #Hbaby = model.components_
    model = NMF(n_components=2, max_iter=10000, tol= 1*10**-8, solver= 'cd', init='nndsvd')#, init='custom')
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IRMix)
    H = model.components_
   # print(W)
    HO = H.copy()
    #print ("H", HO)
    
    H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    #print ("H adjusted", H)
    #print(H)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
    #print(model.fit(IrMix))
    W2 = np.matmul(W,H)
    W2 =np.transpose(W2)
    
   # HO = np.transpose(HO)
   # print(HO)
    #print ("wait", nmfMatcher(IROG,W)[:])
    difsum = 0
# =============================================================================
    for entri in nmfMatcher(IROG,W):
        
        # print('entri', entri)
         
         
        
         #for proport in (H[entri[1][0]]):

        # print("full", (max(HO[entri[0][0]])))
       
         plt.plot(np.linspace(0,1000,numPoints),IROG[entri[0][0],:],color="red")
                           
          #plt.plot(np.linspace(0,1000,numPoints),((W[:,entri[1][0]]*proport)))
         plt.plot(np.linspace(0,1000,numPoints), W[:,entri[1][0]]*np.max(HO))
         plt.gca().invert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()
         DifferenceSum=0
         for n in range (numPoints):
         
             DifferenceSum += abs(IROG[entri[0][0],n]-W[n,entri[1][0]]*np.max(HO))
         difsum += DifferenceSum/numPoints
    print(difsum/2)
         

 
# =============================================================================
nmf2TesterMix(1000,2)




    

# firstRun = pd.DataFrame(data={})
# for n in range (100):
#     firstRun= firstRun.append(NMF2TestRapid(1000))
# print(firstRun)
# plt.scatter(firstRun[0], firstRun[1])
# plt.xlabel("fraction")
# plt.ylabel("error")
# plt.title("Fraction vs Error")
#firstRun.to_csv('SecondRun.csv')

#impure components
#pure components A
#pure B

#forget pure components off
#how many mixtures vs error
#nnvda
#same size

'''
         fractBoi = 0
         lenBoi = 0
         for n in range(len(IROG[entri[0][0],:])):
             
             if (IROG[entri[0][0],n] != 0) and (W[n,entri[1][0]] != 0):
                 fractBoi += IROG[entri[0][0],n]/ W[n,entri[1][0]]
                 #print (IROG[entri[0][0],n]/ W[n,entri[1][0]])
                 lenBoi +=1.0
         fract = fractBoi/lenBoi
         print("fract", fract)
             
         #frac = (IROG[entri[0][0],:])/(W[:,entri[1][0]]*np.max(HO))
        # print(frac)'''        