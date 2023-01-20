#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 13:47:50 2022

@author: abnerkahan
"""

import matplotlib.pyplot as plt
import numpy as np

Ir = np.zeros((5,5 ))

Ir[2,:] = [2,6,9,3,1]

Ir[4,:] = [2,11,8,8,1]
print(Ir)
mask = Ir[:,0] != [0] *len(Ir[0])
print(mask)
print(Ir[mask])