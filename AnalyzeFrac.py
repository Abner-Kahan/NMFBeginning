#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 18:19:14 2022

@author: abnerkahansmack
"""
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF

BigTable = np.load('BigIterFractions.npy')
#print (np.mean(BigTable, axis= 0))
print(np.mean(BigTable[: , 5:12:2], axis=0))
print(np.mean(BigTable[: , 13:20:2], axis=0))
#print(BigTable[: , 13:20:2])
print("Stop \n\n\n")
print(np.mean (np.mean(np.absolute(BigTable[: , 5:12:2]), axis=0)))
print(np.mean(np.mean(np.absolute(BigTable[: , 13:20:2]), axis=0)))
#print(BigTable[11])
#print (np.mean(BigTable[: , 13:20:2]))