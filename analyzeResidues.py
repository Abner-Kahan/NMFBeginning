#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 13:19:41 2022

@author: abnerkahansmack
"""
import numpy as np
H = np.load('tempCompons.npy')
print(np.argmin(H, 1))
H_sort = np.sort(H, 1)
for val in H_sort[0,:10]:
    print(np.where(H == val))