#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 18:03:38 2022

@author: abnerkahan
"""
import matplotlib.pyplot as plt
fig, axs = plt.subplots(4)
axs[0].plot(range(50))
axs[1].plot(range(30))
axs[2].plot(range(20))
axs[3].plot(range(10))
plt.show()

