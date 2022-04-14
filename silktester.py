#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 17:47:47 2022

@author: abnerkahan
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random

import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

bob = np.zeros((3,4))
bob [0,:] = [10,20,60,10]
bob [1,:] = [5,7,54,2]
bob [2,:] = [21,12,7,8]
bob2 = np.apply_along_axis(lambda l : l/np.sum(l),1,bob)
print(bob)
print(bob2)