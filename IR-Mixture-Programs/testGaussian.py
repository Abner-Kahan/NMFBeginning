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
    return Ir
blob = addIRS(10,10000)
plt.plot(np.linspace(0,1000,10000),blob)
plt.clf()
blob = addIRS(10,10000)
plt.plot(np.linspace(0,1000,10000),blob)
