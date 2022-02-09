import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from tempfile import TemporaryFile

from scipy.signal import find_peaks

PeakList = ['BetaTurn', 'Alpha Helix', 'Random Coil', 'Beta Sheet' ]
h = np.load("H.npy")
hnorm = np.load("Hnorm.npy")
hPanda = pd.DataFrame(columns = [1,2,3,4,5])

Humdities = [5,10,20,30,40,50,60,70,80,90,95]

for n in range (len(h[0])):
     hPanda.loc[n]=  (h[:,n] / np.sum(h[:,n]))

#hPanda.index.set_names([])         
print(hPanda)
plt.plot(Humdities,hPanda[:11])
plt.title("Untreated Spectra")
plt.legend([1,2,3,4,5])
plt.ylim(top=.65)
plt.xlabel('% humidity')
plt.show()

plt.clf()
plt.plot(Humdities, hPanda[11:22])
plt.title("MeOH Spectra")
plt.legend([1,2,3,4,5])
plt.ylim(top=.65)
plt.xlabel('% humidity')
plt.show()
plt.clf()

plt.plot(Humdities,hPanda[22:])
plt.title("WA45 Spectra")
plt.legend([1,2,3,4,5])
plt.ylim(top=.65)
plt.xlabel('% humidity')
plt.show()