#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import cauchy
import random
from scipy import constants
#import radis
ReducedMass = 

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
        Ir[start:end]= -1*cauchy.pdf(np.linspace(-1*thickness2,thickness2,end-start))*peakHeight
    return Ir

#IR2Stack= addIRS(30,20000).append[(addIRS(30,20000))]
IR2Stack = np.column_stack((addIRS(30,20000),addIRS(30,20000)))

IR2Stack = np.column_stack((IR2Stack, IR2Stack[:,1] + IR2Stack[:,0]))


plt.plot(np.linspace(4000,400,20000),IR2Stack[:,1] )
plt.plot(np.linspace(4000,400,20000),IR2Stack[:,0] )

plt.plot(np.linspace(4000,400,20000),IR2Stack[:,2] )
#plt.savefig("test.png")




#old code that may not work
#thickness= 10**(random.randrange(10,50)/10 )
#L.extend(-1*cauchy.pdf(np.linspace(-5,5,200)))
#plt.axis([-10,20,-.34,0])
#plt.show()
'''
Irs=[]

for n in range (PointNum):
    Irs+=[0]
    '''
    #plt.plot(np.linspace(4000,400,20000),IR2Stack[:,0] )

#plt.plot(np.linspace(4000,400,20000),IR2Stack[:,1] )
