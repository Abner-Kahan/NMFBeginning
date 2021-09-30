#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
import scipy.stats as stats
import random
#from scipy import constants
from sklearn.decomposition import NMF
#import seaborn as sns
#import matplotlib.ticker as mtick
def nmfListTester(spectra):
    #IR_stack = addIRS(10,20000)
    #plt.plot(np.linspace(0,1000,20000),IR_stack)
    #plt.gca().invert_yaxis()
    #plt.title("0 of " + str(spectra)+ " : Original Spectra")
    #plt.savefig("0of" + str(spectra)+ ":Original Spectra.png")
   # plt.show()
    IR_stack = np.empty(20000)
    ind = 0
    for n in spectra:
        n = int(n * wholeNum(spectra, 1))

        for q in range (n):
            IR_stack = np.column_stack((IR_stack,addIRS(10,20000)))
        
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,ind])
        plt.gca().invert_yaxis()
        plt.title(str(n+1)+ " of " + str(spectra)+ ": Original Spectra")
        #plt.savefig((str(n+1)+ "of" + str(spectra)+ "_ Original Spectra.png"))
        plt.show()
        plt.close()
        ind +=1
    plt.plot(np.linspace(0,1000,20000),IR_stack)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Original Spectra")
    plt.xlabel("cm^-1")
   # plt.savefig((str(spectra)+ " Original Spectra.png"))
    plt.show()
    model = NMF(n_components=len(spectra), init='random', random_state=0, shuffle=1 )
    
    W = model.fit_transform(IR_stack)
    for n in range(len(spectra)):
        plt.plot(np.linspace(0,1000,20000),W[:,n],markersize=1)
        plt.gca().invert_yaxis()
        plt.title(str(n)+ " of " + str(spectra)+ ": Calculated Spectra")
       # plt.savefig((str(n)+ " of " + str(spectra)+ "_ Calculated Spectra.png"))
        plt.show()        
        plt.close()
    plt.plot(np.linspace(0,1000,20000),W)
    plt.gca().invert_yaxis()
    plt.title(str(spectra)+ " Calculated Spectra")
    plt.xlabel("cm^-1")
    plt.close()
    for entri in nmfMatcher(IR_stack,W):
        plt.plot(np.linspace(0,1000,20000),IR_stack[:,entri[0][0]],color="red")
        plt.plot(np.linspace(0,1000,20000),W[:,entri[1][0]])
        plt.legend(["Original", "Calculated"])
        plt.gca().invert_yaxis()
        plt.title(str(entri[0][0])+ " Both Spectra")
        plt.xlabel("cm^-1")
        plt.show()
        plt.close()