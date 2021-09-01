#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import cauchy
xbro= np.linspace(-.99,.99,10000),
#randList= [200,100]
plt.plot(xbro, cauchy.pdf(xbro))
#cauchy.pdf
#plt.plot(np.linspace(cauchy.pdf(0.01),
               # cauchy.pdf(0.99), 100))