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
    
    
def nmf2TesterMix(numPoints):
    fraction1 = random.random()
    print(f'The expected fraction is {fraction1:.3}')
    IRF = np.zeros((2,numPoints))
    IR0 = addIRS(8,numPoints)
    IrPlotter(IR0,"FirstPlot")
    IR1 = addIRS(8,numPoints)
    IrPlotter(IR1,"SecondPlot")
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((3,numPoints))
    IrMix[0,:]=IR0
    IrMix[1,:]=IR1
    IrMix[2,:] = IR0*fraction1 + IR1*(1-fraction1) 
    #IrMix[3,:] = IR0*(1-fraction2)  + IR1*fraction1
    IrPlotter( IrMix[0,:],"FirstMix")
    IrPlotter(IrMix[1,:],"SecondMix")
    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd',  max_iter=15000, tol= 1*10**-7)
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    H = model.components_
    print(W)
    HO = H.copy()
    print(H)
    H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    #print(H)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
    #print(model.fit(IrMix))
# =============================================================================
    for entri in nmfMatcher(IRF,W):
        
         plt.plot(np.linspace(0,1000,numPoints),IRF[:,entri[0][0]],color="red")
         if H[0,0]>.01:
             print(f'The calculated fraction of the first is {H[0,2]:.3}.')
             print(f'The calculated fraction of the second is {H[1,2]:.3}.')

             
         else:
             print(f'The calculated fraction of the first is {H[1,2]:.3}')
             print(f'The calculated fraction of the second is, {H[0,2]:.3}.')
         plt.plot(np.linspace(0,1000,numPoints),(W[:,entri[1][0]]*(max(HO[entri[1][0]]))))
        # print("full", (max(HO[entri[0][0]])))
        


         plt.gca().invert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()
# =============================================================================

def nmfTesterMix(fracList,numPoints):
    IRF = np.zeros((2,numPoints))
    IR0 = addIRS(8,numPoints)
    IrPlotter(IR0,"Plot A")
    IR1 = addIRS(8,numPoints)
    IrPlotter(IR1,"Plot B")
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((len(fracList),numPoints))
    ind = 0
    for frac in fracList:
        IrMix[ind] = IR0*frac + IR1*(1-frac) 
       #IrPlotter( IrMix[ind,:], "Mix Number "+ str(ind)+ "  "+ str(frac*100)+ "% A "+ str((1-frac)*100) + "% B")
        IrPlotter( IrMix[ind,:], f'Mix Number {ind} {frac*100:,.1f}% A {(1-frac)*100:,.1f}% B')
        ind+=1
        
       

    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd',  max_iter=15000, tol= 1*10**-7)
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    H = model.components_
    HO = H.copy()
    #print(H)
    H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    strin = 'The percentages of A are: '
    for n in H[0]:
        strin += str(np.round(n*100,1))
        strin +='%, '
    print(strin) 
    strin2 = 'The percentages of B are: '
    for n in H[1]:
        strin2 += str(np.round(n*100,1))
        strin2 +='%, '
    print(strin2)    
    
            
            
      #  print 
    #print("real one")
    #print(f'The percentages of A are {H[0]:.2%f}')
    #print(H)
    IrPlotter(W[:,0], "First Calc Spectra")
    IrPlotter(W[:,1], "Second Calc Spectra")
    #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
    #print(model.fit(IrMix))
# =============================================================================
    W = np.matmul(W,H)
    for entri in nmfMatcher(IRF,W):
        
         plt.plot(np.linspace(0,1000,numPoints),IRF[:,entri[0][0]],color="red")
        # if H[0,0]>.01:
            # print("The fraction of the first is", H[0,2])
             #print("The fraction of the second is", H[1,2])

             
        # else:
             #print("The fraction of the first is", H[1,2])
             #print("The fraction of the second is", H[0,2])
         plt.plot(np.linspace(0,1000,numPoints),(W[:,entri[1][0]]*(max(HO[entri[1][0]]))))
        # print("full", (max(HO[entri[0][0]])))
        


         plt.gca().inNMF22vert_yaxis()
         plt.legend(["Original", "Calculated"])
         plt.title(str(entri[0][0])+ " Both Spectra")
         plt.xlabel("cm^-1")
         plt.show()
         plt.clf()

#nmfTesterMix([1,.5,.25,0,.8,.3],10000)
# =============================================================================
#nmf2TesterMix(10000)

def NMF2TestRapid(numPoints):
    fraction1 = random.random()
    
    IRF = np.zeros((2,numPoints))
    IR0 = addIRS(8,numPoints)
    IR1 = addIRS(8,numPoints)
    IRF[0,:] = IR0
    IRF[1,:] = IR1
    IRF= np.transpose(IRF)
    
    IrMix = np.zeros((3,numPoints))
    IrMix[0,:]=IR0
    IrMix[1,:]=IR1
    IrMix[2,:] = IR0*fraction1 + IR1*(1-fraction1) 
    #IrMix[3,:] = IR0*(1-fraction2)  + IR1*fraction1
    IrMix= np.transpose(IrMix)
    model = NMF(n_components=2, init='nndsvd',  max_iter=15000, tol= 1*10**-7)
    #it seems that mu gives more close results
    #must analyze errors and create plots
    W = model.fit_transform(IrMix)
    H = model.components_
    W2 = np.matmul(W,H)
    error = np.sum(abs(IrMix-W2))/IrMix.size
    if error > 1*10**-7:
        HO = H.copy()
        H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
        for entri in nmfMatcher(IRF,W):
            plt.plot(np.linspace(0,1000,numPoints),IRF[:,entri[0][0]],color="red")
            if H[0,0]>.01:
                print(f'The calculated fraction of the first is {H[0,2]:.3}.')
                print(f'The calculated fraction of the second is {H[1,2]:.3}.')
            else:
              print(f'The calculated fraction of the first is {H[1,2]:.3}')
              print(f'The calculated fraction of the second is, {H[0,2]:.3}.')
            plt.plot(np.linspace(0,1000,numPoints),(W[:,entri[1][0]]*(max(HO[entri[1][0]]))))
            print("full", (max(HO[entri[0][0]])))
        plt.gca().invert_yaxis()

        plt.legend(["Original", "Calculated"])
        plt.title(str(entri[0][0])+ " Both Spectra")
        plt.xlabel("cm^-1")
        plt.show()
        plt.clf()
    return pd.DataFrame(data=[[fraction1, error]])
    #print(W)
    #HO = H.copy()
    #print(H)
   # H = np.apply_along_axis(lambda l :l/np.amax(l) ,1,H)
    
#     #print ("taco", IrMix[:,0:10])
#     print ("squash", W2[:,0:10])
#     ##print(H)
# #   IrPlotter(W[:,0], "First Calc Spectra")
#    # IrPlotter(W[:,1], "Second Calc Spectra")
#     #print (np.mean(np.where(W[:,1]>0))/np.mean((np.where(W[:,0]>0)))
#     #print(model.fit(IrMix))
# # =============================================================================
#     for entri in nmfMatcher(IRF,W):
        
#           #plt.plot(np.linspace(0,1000,numPoints),IRF[:,entri[0][0]],color="red")
#           if H[0,0]>.01:
#               print(f'The calculated fraction of the first is {H[0,2]:.3}.')
#               print(f'The calculated fraction of the second is {H[1,2]:.3}.')

             
#           else:
#               print(f'The calculated fraction of the first is {H[1,2]:.3}')
#               print(f'The calculated fraction of the second is, {H[0,2]:.3}.')
#          #plt.plot(np.linspace(0,1000,numPoints),(W[:,entri[1][0]]*(max(HO[entri[1][0]]))))
#         # print("full", (max(HO[entri[0][0]])))
        


    

firstRun = pd.DataFrame(data={})
for n in range (100):
    firstRun= firstRun.append(NMF2TestRapid(10000))
print(firstRun)
plt.scatter(firstRun[0], firstRun[1])
firstRun.to_csv('SecondRun.csv')
