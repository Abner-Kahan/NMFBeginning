#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 17:17:27 2022

@author: abnerkahan
"""

paramlist = ['random-cd-frobenius' ,
               'nndsvd-cd-frobenius' ,
               'nndsvda-cd-frobenius' ,
               'nndsvdar-cd-frobenius' , 
               'random-mu-frobenius' ,
               'nndsvd-mu-frobenius' ,
               'nndsvda-mu-frobenius' ,
               'nndsvdar-mu-frobenius' ,
               'random-mu-kullbackLeibler' ,
               'nndsvd-mu-kullbackLeibler' ,
               'nndsvda-mu-kullbackLeibler' ,
               'nndsvdar-mu-kullbackLeibler'  ]
for par in paramlist:
    print(par[:par.find('-')])
    print(par[par.find('-')+1:par.rfind('-')])
    print(par[par.rfind('-')+1:])
    