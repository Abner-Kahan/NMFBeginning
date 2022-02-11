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
import subprocess
from scipy.signal import find_peaks




def fetchIr(path):
    logfile =  open(path, 'r')
    lines = logfile.readlines()
    logfile.close
    lines = filter(lambda x:'#' in x and '.' in x, lines)
    lines = list(lines)
    spectra = np.zeros((len(lines),2))
    index = 0
    for line in lines:
        spectra[index] = line.split()[1:3]
        index+=1
                                      
   
    return spectra

#print(fetchIr('3A6/gas/input_ir.txt'))
#print(fetchIr('2A6/gas/input_ir.txt'))



def IrPlotter(item,title,ran1,ran2,leg = [], multiple = False):
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,len(item)),item,markersize=.1)
    else:
        for n in item:
            plt.plot(np.linspace(ran1,ran2,len(n)),n,markersize=.1)
    if len (leg) > 0:
        plt.legend(leg)
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.ylim(0,40000)
    plt.xlabel("cm^-1")
    plt.show()
    plt.clf()    


def gaussian_broadening(spectra, broaden, ran1,ran2,resolution=1):
 
    """ Performs gaussian broadening on IR spectrum
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum
    
    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """
    IR = np.zeros(resolution*(int(ran2-ran1) + 1))
    X = np.linspace(ran1,ran2, resolution*(int(ran2-ran1)+1))

    #for f, i in zip(spectra[:,0]):
      #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
      #  IR=np.vstack((X, IR)).T #tspec
   
    freq = spectra[:,0]
    inten = spectra[:,1]
    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
    return IR
output = subprocess.check_output(('find . -name input_ir.txt'), shell=True)

output = output.split()
print (output)
ind = 0
for n in output:
    n= n.decode('UTF-8')
    output[ind] = n
    ind +=1
for Type in ['/gas/', '/wat_gas/', '/pcm/', '/wat_pcm/' ]:
    type_filter = filter(lambda x: ('helix' in x) and (Type in x), output)
    Inter =list(type_filter)
    specs = []
    legend = []
    for spec in Inter:
       i = spec.index('/',9)
       legend.append(spec[8:i])
       specs.append(gaussian_broadening(fetchIr(spec), 20, 1450, 1700))
       #print("check" ,specs)
    IrPlotter(specs, "Helix " + Type[1:-1] + ' amide region', 1450, 1750, leg = legend, multiple = True)

for Type in ['/gas/', '/wat_gas/', '/pcm/', '/wat_pcm/' ]:
    type_filter = filter(lambda x: ('A6' in x) and (Type in x), output)
    Inter =list(type_filter)
    specs = []
    legend = []
    for spec in Inter:
       #i = spec.index('/',9)
       legend.append(spec[2:5])
       specs.append(gaussian_broadening(fetchIr(spec), 20, 1450, 1700))
       #print("check" ,specs)
    IrPlotter(specs, "Bsheets " + Type[1:-1] + ' amide region', 1450, 1750, leg = legend, multiple = True)

    
   # 
    
# =============================================================================
# 
# IrPlotter( gaussian_broadening(fetchIr('2A6/gas/input_ir.txt'), 20, 0, 4000),"2A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('2A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('3A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('3A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('3A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('3A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('3A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# IrPlotter( gaussian_broadening(fetchIr('3A6/gas/input_ir.txt'), 20, 0, 4000),"3A6 gas", 0, 4000)
# =============================================================================
