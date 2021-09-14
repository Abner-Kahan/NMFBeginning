#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import cauchy
import random


Irs=[]
PointNum=20000
for n in range (PointNum):
    Irs+=[0]

def addIRS(peakNum,L):
    for n in range(peakNum):
        xloc=  random.randrange(len(L))
        #thickness= 10**(random.randrange(10,50)/10 )
        thickness = random.randrange(10,1000)
        thickness2 = random.randrange(1,5)
        peakHeight=random.random()
#add width
        if (xloc + thickness+1) > PointNum:
            end = PointNum
        else: 
            end = (xloc + thickness+1)
        if (xloc-thickness < 0):
            start= 0
        else:
            start = xloc - thickness
        L[start:end]= -1*cauchy.pdf(np.linspace(-1*thickness2,thickness2,end-start))*peakHeight
    return L
addIRS(20,Irs)

#L.extend(-1*cauchy.pdf(np.linspace(-5,5,200)))
plt.plot(np.linspace(4000,400,PointNum),Irs)
plt.savefig("test.png")

#plt.axis([-10,20,-.34,0])
#plt.show()


