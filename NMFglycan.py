# -*- coding: utf-8 -*-
"""
Created on Tue May 24 09:57:49 2022

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

mapGL2 =  open('100inc.txt', 'r')
maptext = mapGL2.read()  
print(maptext)
mapGL2.close()
#print(maptext)
print(type(maptext))
#print(type(maptext[0]))
contactList = maptext.split("{")
print(len(contactList))
