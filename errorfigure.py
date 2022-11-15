#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 15:58:19 2022

@author: abnerkahan
"""
import random
def rotate(l, n):
    return l[n:] + l[:n]
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.spatial import distance
fig, axs = plt.subplots(2)
P =[]
for n in range(8):
    P.append(random.random()/2)
#P = [.05, .1, .2, .05, .15, .25, .08, .12]
Q =[]
for n in P:
  Q.append(n+.1)
KLError = sum(scipy.special.kl_div(P,Q))
EucError = distance.euclidean(P,Q)
axs[0].plot(P, color = 'orange')
axs[0].plot(Q, color = 'blue')

axs[0].set_title(f"Frobenius Error {EucError:.3} > KL Error {KLError:.3}",
                 
                 fontfamily='serif', loc='left', fontsize='medium')
plt.subplots_adjust(hspace=(.4) )        
axs[0].legend(['Original', 'Y-shifted upwards'])
P2 = [.05, .1, .2, .05, .15, .25, .08, .12]
#Q2 =rotate(P2,4)

Q2 = rotate(P,4)


KLError2 = sum(scipy.special.kl_div(P,Q2))
EucError2 = distance.euclidean(P,Q2)
axs[1].plot(P, color = 'orange')
axs[1].plot(Q2, color = 'green')
axs[1].set_title(f"KL Error {KLError2:.3} > Frobenius Error {EucError2:.3}",fontfamily='serif', loc='left', fontsize='medium')
axs[1].legend(['Original', 'Original Shifted by 4'])

