# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 10:20:57 2022

@author: Abner
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import savgol_filter
from scipy.signal import find_peaks


amide1_ran1 =1600
amide1_ran2 =1715


Humdities = [5,10,20,30,40,50,60,70,80,90,95]
solvents = ["Untreated", "MeOH", "WA45"]
broad = 15
x_range = np.linspace(amide1_ran1,amide1_ran2,(amide1_ran2-amide1_ran1+1))

def fetchIr(path):
    logfile =  open(path, 'r')
    logtest = logfile.readlines()
   
    logfile.close()
    
    rows = len(logtest)
    columns = 2
              
    IrArray = np.zeros((rows,columns ))
    
    x = 0
    for line in logtest:
        line2 = line.split()
        IrArray[x] = line2
        x +=1

    
    
    
    return IrArray


#fetch
#print(fetchIr('UntreatedSample.txt',3))
print(fetchIr('BGH2c_h.dat'))
def IrPlotter(item,title,ran1,ran2, leg = [], multiple = False):
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,len(item)),item,markersize=.1)
    else:
        for n in item:
            plt.plot(np.linspace(ran1,ran2,len(n)),n,markersize=.1)
    if len (leg) > 0:
        plt.legend(leg,fontsize='x-small')
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")    
    plt.show()
    plt.clf() 
    
    
plt.plot(fetchIr('Fuc_a16a.dat')[:,0],fetchIr('Fuc_a16a.dat')[:,1])
plt.show()
plt.clf()
plt.plot(fetchIr('Fuc_a16b.dat')[:,0],fetchIr('Fuc_a16b.dat')[:,1])
plt.show()
plt.clf()

plt.plot(fetchIr('Fuc_aLex.dat')[:,0],fetchIr('Fuc_a16a.dat')[:,1])
plt.show()
plt.clf()
plt.plot(fetchIr('Fuc_bLex.dat')[:,0],fetchIr('Fuc_a16b.dat')[:,1])
plt.show()
plt.clf()

plt.plot(fetchIr('Fuc_aBGH2.dat')[:,0],fetchIr('Fuc_a16a.dat')[:,1])
plt.show()
plt.clf()
plt.plot(fetchIr('Fuc_bBGH2.dat')[:,0],fetchIr('Fuc_a16b.dat')[:,1])
plt.show()
plt.clf()
