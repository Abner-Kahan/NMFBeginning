#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import cauchy
import random


Irs=[]
for n in range (6000):
    Irs+=[0]

def addIRS(peakNum,L):
    for n in range(peakNum):
        xloc=  random.randrange(len(L))
#add width
        if (xloc + 201) > 6000:
            end = 6000
        else: 
            end = (xloc + 201)
        if (xloc-200 < 0):
            start= 0
        else:
            start = xloc - 200
        L[start:end]= -1*cauchy.pdf(np.linspace(-5,5,end-start))
    return L
addIRS(10,Irs)

#L.extend(-1*cauchy.pdf(np.linspace(-5,5,200)))
plt.plot(np.linspace(0,1000,6000),Irs)
#plt.axis([-10,20,-.34,0])
#plt.show()


