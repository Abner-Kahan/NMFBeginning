# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 16:32:22 2022

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

Humdities = [5,10,20,30,40,50,60,70,80,90,95]

H_n = np.load('Hnorm.npy')
H = np.load('H.npy')
plt.plot(Humdities, H[0, :11], label = "1620,1689")
plt.plot(Humdities, H[1, :11], label = "1662")
plt.plot(Humdities, H[2, :11], label = "1639")
plt.plot(Humdities, H[3, :11], label = "1644")
plt.legend()