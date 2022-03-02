#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 18:03:45 2022

@author: abnerkahansmack
"""
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import random
def coolF(x,k):
    return k*x**2
def superCoolF(x,m,b):
    return m*x+b
ySolvedata = []
ySolvedata2 = []
ydata = [10,35,93,169,283,347,490,629,800]

fit = curve_fit(coolF,list(range(1,10)) , ydata)
fit2 = curve_fit(superCoolF,list(range(1,10)) , ydata)
print(fit2)
for n in range(1,10):
    ySolvedata +=[coolF(n, fit[0])]
    ySolvedata2 +=[superCoolF(n, *fit2[0])]
    
plt.plot(range(1,10), [10,35,93,169,283,347,490,629,800])
plt.plot(list(range(1,10)),ySolvedata, color='red')
plt.plot(list(range(1,10)),ySolvedata2, color='black')
print(fit)

