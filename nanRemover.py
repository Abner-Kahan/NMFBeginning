#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:33:36 2022

@author: abnerkahan
"""
import numpy as np
data = np.array([[4, np.nan, 6, np.nan, 10, 11, 14, 19, 22], [np.nan, 2, 6, 10, 10, 11, 14, 19, 22]])
print(data)
print(data[~np.isnan(data)])