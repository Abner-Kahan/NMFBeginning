#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate as integrate
import scipy.special as special
import random
#from scipy import constants
from sklearn.decomposition import NMF

#np.set_printoptions(precision=3)

def IrPlotter(item,title,colori):
    plt.plot(np.linspace(0,1000,len(item)),item,markersize=.1,color=colori)
    plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    

   
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
Spectra1=addIRS(10,4000)
Spectra2=addIRS(10,4000)

 
DifferenceSum=0
for n in range (4000):
    DifferenceSum += abs(Spectra1[n]-Spectra2[n])
print(DifferenceSum/4000)
    


IrPlotter(Spectra1,"First Spectra",'b')
IrPlotter(Spectra2,"Second Spectra",'y')
IrPlotter(abs(Spectra1-Spectra2),"Difference","magenta")
plt.legend()

plt.show()


