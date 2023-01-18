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
    #lines =sorted(lines)
    spectra = np.zeros((len(lines),2))
    index = 0
    for line in lines:
        spectra[index] = line.split()[1:3]
        index+=1
                                      
    #print(spectra)
    return spectra

#print(fetchIr('3A6/gas/input_ir.txt'))
#print(fetchIr('2A6/gas/input_ir.txt'))



def IrPlotter(item,title,ran1,ran2,leg = [], multiple = False):
    #colors = np.array([(254,230,206), (253,174,107),(230,85,13)])
    #colors = colors/256
    colors =['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']
    print(colors)

        
    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),item,markersize=.1)
    else:
        index = 0
        for n in item:
            plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),n,markersize=.1,color=colors[index+1])
            index += 1
    if len (leg) > 0:
        plt.legend(leg)
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.ylim(0,45000)
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
   
    freq = spectra[:,0]*.965
   # print(np.max(freq))
    #print(np.min(freq))
    inten = spectra[:,1]
    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
        
    
    return IR
OperatingSystem= 'Windows'
ind = 0


   # output = subprocess.check_output(('dir /S  /B input_ir.txt'), shell=True)
   # output = output.replace(b'\\',b'/')
   # output = output.replace(b'\n',b'')
    #output = output.replace(b'\r',b'')
    
    #root = subprocess.check_output('PWD', shell=True)
    #print(root)
    #output = output.split(B'C:/Users/Abner/Documents/Python Scripts')

    #output=output[::2]
output = glob.glob('*A6/**/input_ir.txt',recursive = True)
outputWin = glob.glob('*A6\\**\\input_ir.txt',recursive = True)
print (output)

#print (output) 
if OperatingSystem =='Windows':
    TypeList = ['\\gas\\', '\\wat_gas\\', '\\pcm\\', '\\wat_pcm\\']
    for Type in TypeList:
        type_filter = filter(lambda x: ('A15' in x) and (Type in x), output)
        Inter =list(type_filter)
        
    #print(Inter)
        specs = []
        legend = []
    for spec in Inter:
       i = spec.index('\\',9)
       legend.append(spec[8:i])
       specs.append(gaussian_broadening(fetchIr(spec), 20, 1450, 1800))
       #print("check" ,specs)
       #print('specs', len(specs))
    #print(specs)
       
    IrPlotter(specs, "Helix " + Type[1:-1] + ' amide region', 1450, 1800, leg = legend, multiple = True)

    for Type in TypeList:
        type_filter = filter(lambda x: ('A15' in x) and (Type in x), output)
        Inter =list(type_filter)
        specs = []
        legend = []
        for spec in Inter:
           i = spec.index('/',10)
           legend.append(spec[:i])
           specs.append(gaussian_broadening(fetchIr(spec), 20, 1450, 1800))
           #print("check" ,specs)
        IrPlotter(specs, "helix" + Type[1:-1] + ' amide region', 1450, 1800, leg = legend, multiple = True)
    else:
        TypeList = ['/gas/', '/wat_gas/', '/pcm/', '/wat_pcm/']
        
     
if OperatingSystem =='Linux':
    TypeList = ['/gas/', '/wat_gas/', '/pcm/', '/wat_pcm/']
# =============================================================================
#     for Type in TypeList:
#         type_filter = filter(lambda x: ('A15' in x) and (Type in x), output)
#         Inter =list(type_filter)
#         print("test", Inter)
#         specs = []
#         legend = []
# 
#         
#         for spec in Inter:
#            #m = spec.index('/',4) 
#            #i = spec.index('/',9)
#            legend.append(spec[:])
#        #i = spec.index('/',9)
#        #legend.append(spec[8:i])
#            specs.append(gaussian_broadening(fetchIr(spec), 20, 1450, 1800))
#        #print("check" ,specs)
#        #print('specs', len(specs))
#     #print(specs)
#         IrPlotter(specs, "Helix " + Type[1:-1] + ' amide region', 1450, 1800, leg =legend, multiple = True)
# =============================================================================

    for Type in TypeList:
        type_filter = filter(lambda x: ('A6' in x) and (Type in x), output)
        Inter =list(type_filter)
        Inter =sorted(Inter)
        specs = []
        legend = []
        for spec in Inter:
           #m = spec.index('/',4) 
           #i = spec.index('/',8)
           legend.append(spec[:])
           specs.append(gaussian_broadening(fetchIr(spec), 20, 1450, 1800))
           #print("check" ,specs)
        IrPlotter(specs, "Bsheets " + Type[1:-1] + ' amide region', 1450, 1800, leg = legend, multiple = True)
    else:
        TypeList = ['/gas/', '/wat_gas/', '/pcm/', '/wat_pcm/']        




    
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
