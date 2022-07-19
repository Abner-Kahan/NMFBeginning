#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:19:41 2022

@author: abnerkahansmack
"""
import numpy as np
import matplotlib.pyplot as plt

H = np.load('tempCompons.npy')
print(np.argmin(H, 1))
H_sort = np.sort(H, 1)
for val in H_sort[0,:10]:
    print(np.where(H == val))
for search in range (0,10000,10):
#search= 
    print(H[0, search], search)
#print(search)
plt.plot(np.linspace(0, len(H[0,:]), len(H[0,:])), H[0,:], linewidth=.3)