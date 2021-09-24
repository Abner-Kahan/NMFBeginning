#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from scipy.stats import cauchy
import random
#from scipy import constants
from sklearn.decomposition import NMF


def addIRS(peakNum,PointNum):
    Ir = np.zeros(PointNum)
    for n in range(peakNum):
        xloc=  random.randrange(PointNum)
        
        thickness = random.randrange(10,1000)
        thickness2 = random.randrange(1,5)
        peakHeight=random.random()
#This is supposed to deal with peaks toward the ends of the spectra
        if (xloc + thickness+1) > PointNum:
            end = PointNum
        else: 
            end = (xloc + thickness+1)
        if (xloc-thickness < 0):
            start= 0
        else:
            start = xloc - thickness
        Ir[start:end]= cauchy.pdf(np.linspace(-1*thickness2,thickness2,end-start))*peakHeight
    return Ir

#IR2Stack= addIRS(30,20000).append[(addIRS(30,20000))]
IR2Stack = np.column_stack((addIRS(30,20000),addIRS(30,20000)))
#plt.plot(np.linspace(4000,400,20000),IR2Stack)
#plt.title("Two Original Spectra")
#plt.xlabel("cm^-1")
#plt.savefig("basespectra.png")
#plt.show()
#plt.close()
#IR2Stack = np.column_stack((IR2Stack, IR2Stack[:,1] + IR2Stack[:,0]))


#plt.plot(np.linspace(4000,400,20000),IR2Stack[:,1] )
#plt.plot(np.linspace(4000,400,20000),IR2Stack[:,0] )

#plt.plot(np.linspace(4000,400,20000),IR2Stack[:,2] )
#plt.savefig("test.png")
'''
model = NMF(n_components=2, init='nndsvda', random_state=0)
W = model.fit_transform(IR2Stack)
plt.plot(np.linspace(4000,400,20000),W[:,0],color = 'cyan')
plt.show()
plt.close()
plt.plot(np.linspace(4000,400,20000),W[:,1],color = 'yellow')
plt.title("Two Calculated Spectra")
plt.xlabel("cm^-1")
#plt.savefig("calcSpectra.png")
plt.show()
plt.close()
H = model.components_
'''
def nmfTester(spectra):
    IR_stack = addIRS(30,20000)
    plt.plot(np.linspace(4000,400,20000),IR_stack)
    plt.title("0 of " + str(spectra)+ " : Original Spectra")
    plt.savefig("0of" + str(spectra)+ ":Original Spectra.png")
    plt.show()
    plt.close()
    for n in range(spectra - 1):
        IR_stack = np.column_stack((IR_stack,addIRS(30,20000)))
        plt.plot(np.linspace(4000,400,20000),IR_stack[:,n+1])
        plt.title(str(n+1)+ " of " + str(spectra)+ ": Original Spectra")
        plt.savefig((str(n+1)+ " of " + str(spectra)+ "_ Original Spectra.png"))
        plt.show()
        plt.close()
    plt.plot(np.linspace(4000,400,20000),IR_stack)
    plt.title(str(spectra)+ " Original Spectra")
    plt.xlabel("cm^-1")
    plt.savefig((str(spectra)+ " Original Spectra.png"))
    plt.show()
    model = NMF(n_components=spectra, init='random', random_state=0,alpha=1,max_iter=3000 )
    W = model.fit_transform(IR2Stack)
    for n in range(spectra):
        plt.plot(np.linspace(4000,400,20000),W[:,n],markersize=1)
        plt.title(str(n)+ " of " + str(spectra)+ ": Calculated Spectra")
        plt.savefig((str(n)+ " of " + str(spectra)+ "_ Calculated Spectra.png"))
        plt.show()
        plt.close()
    plt.plot(np.linspace(4000,400,20000),W)
    plt.title(str(spectra)+ " Calculated Spectra")
    plt.xlabel("cm^-1")
    return W
nmfTester(1)
        

