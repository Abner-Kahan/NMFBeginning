#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 13:37:57 2021

@author: abnerkahan
"""
import matplotlib.pyplot as plt
import numpy as np
import random
'''

#print(np.random.normal(500,10,100))
#plt.plot(np.linspace(0,1000,10000),np.random.normal(500,10,10))
Ir = np.zeros(1000)
#print(random.choice(np.where(Ir == 0)[0]))
Ir[100:900] = 1
print(random.choice(np.where(Ir == 0)[0]))

Irlocs = np.zeros(100)
for b in range(100):
    print(b)
  
p=2
def g():
    try:
        print (1/0)
    except:
        print("no good")
g()
'''
Irlocs = np.zeros(100)
Irlocs [65]= 1
for b in range(80):
    if Irlocs[b:b+5] >0:
        print(b)
        