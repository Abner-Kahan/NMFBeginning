#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import re


logfile =  open('Tri_A131/Tri_A131_0000/input.log', 'r')
logtest = logfile.read()
freqstri =re.findall('Frequencies\D*(\d+.\d+)\D*(\d+.\d+)\D*(\d+.\d+)',logtest)
IRIntenstri =re.findall('IR Inten\D*(\d+.\d+)\D*(\d+.\d+)\D*(\d+.\d+)',logtest)


IrDict =[]
for freqTuple,intTuple in zip(freqstri,IRIntenstri):
    for n,p in zip(freqTuple,intTuple):
        IrDict.append( [float(n), float(p)])

Irs = np.array(IrDict)
Irs[:,1] = 100*Irs[:,1]/np.amax(Irs[:,1])
plt.vlines(Irs[:,0],0,Irs[:,1])
plt.xlabel("Wavenumber cm^-1")
plt.ylabel("Intensity %")
plt.gca().invert_yaxis() 
plt.gca().invert_xaxis()  #print(logtest)
#print(intensities.findall(logtest))
