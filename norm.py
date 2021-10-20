#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

"""
Created on Wed Oct 20 14:18:45 2021

@author: abnerkahan
"""
np.set_printoptions(precision=3)

def normilize(L):
    ma = 0
    for subL in L:
        ma = max(subL+[ma])
        n = 0
        while n < len(subL):
            subL[n] =  subL[n] /  ma
            n += 1
        ma = 0
    return L



print(normilize([[2.92083572, (6.00226641*10**-9), 1.46041786], [0.00000000, 2.54629877, 1.27314938]]))
bob = np.arange(9).reshape(3,3)
print(bob)
def divMin(l):
    l = l/np.amax(l)
    return l
#for n in range (len(bob)):
bob = np.apply_along_axis(lambda l :l/np.amax(l) ,1,bob)
print(bob)

ind = 3
frac =.48
print(f'Mix Number {ind} {frac*100}% A {(1-frac)*100}% B')