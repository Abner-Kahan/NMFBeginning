#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from scipy.stats import cauchy
import random
#from scipy import constants
from sklearn.decomposition import NMF

data = np.array([[1,8,3,3,4],
                 [1,8,9,9,4],
                 [1,8,3,3,4]])
data = np.unique(data,axis=0)
#print(data)


data2= [[1,8,3,3,4],
                 [1,8,9,9,4],
                 [1,8,9,9,4],
                 [1,8,3,3,4]]

for match in data2:
      if data2.count(match) > 1:
          data2.remove(match)
           
print (data2)         