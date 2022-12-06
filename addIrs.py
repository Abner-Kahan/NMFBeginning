#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 12:37:58 2022

@author: abnerkahan
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

def addIRS(PointNum):
    IR = np.zeros((PointNum,4))
    X =  np.linspace(0, PointNum, PointNum)

    thickness =  PointNum * .1

    xloc = [20,45,55,80]
    for n in range(4):
        IR[:,n] = .5*np.exp(-0.5*((X-xloc[n])/int(thickness))**2)



    return IR 

IRS = addIRS(120)
plt.plot(IRS)