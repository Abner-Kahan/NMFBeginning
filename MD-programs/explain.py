#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 10:32:11 2022

@author: abnerkahansmack
"""
import pandas as pd
import numpy as np
tom = pd.DataFrame(np.zeros((9,9)))
tom.columns = ['Residue 1', 'Residue 2','Residue 3','Residue 4','Residue 5','Residue 6','Residue 7','Residue 8','Residue 9' ]
tom.index =[ 'Residue 1', 'Residue 2','Residue 3',,'Residue 4','Residue 5','Residue 6','Residue 7','Residue 8','Residue 9' ]
for n in range(1, 82):
    tom[n] = n
print(n)
