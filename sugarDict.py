#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 11:00:09 2022

@author: abnerkahan
"""


sugarList = ['m5-45d.txt','m6-45d.txt','m7-45d.txt', 'm8-45d.txt', 'm9-45d.txt','n1fb-45d.txt' , 
             'n2-45d.txt', 'n2b-45d.txt','n2f-45d.txt','n33-45d.txt','n36-45d.txt', 's1-45d.txt', 
             's2-45d.txt',  's33-45d.txt', 's36-45d.txt']
residueCount = [7,8, 9 ,10 ,11, 8, 7, 8, 8, 8, 8, 10, 11, 14, 14   ]

sugarDirectory = {}
for sugar in zip(sugarList,residueCount) :
    sugarDirectory[sugar[0][:sugar[0].find('-')]] = sugar[1]
print(sugarDirectory)