#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:16:50 2022

@author: abnerkahan
"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scipy
import scipy.stats as stats
import random
import pdb
import re
import glob
from sklearn.decomposition import NMF
from scipy.signal import savgol_filter
from scipy.signal import find_peaks,deconvolve, peak_widths
from scipy.optimize import curve_fit
import matplotlib.image as mpimg


ran1 =1600
ran2 = 1715
broad = 15
x_range = np.linspace(ran1,ran2,ran2-ran1+1)


def ImaginaryEnergy(spectra):
    peaks, _ = find_peaks(spectra)
    results_half = peak_widths(spectra, peaks, rel_height=0.5)
    ImaginaryEnergy = np.average (results_half[0])/2
    return ImaginaryEnergy
def Y(wavelength, spectra,V_oi):
    L = scipy.misc.derivative (lambda x:(spectra[int(x)]), wavelength)/spectra[int(wavelength)]
    return (L**-1)/(L**-2 + V_oi**2)

def num(wavelength, TheoSpec,ExperSpec,V_oi):
    return (Y(wavelength, TheoSpec,V_oi) - Y(wavelength, ExperSpec,V_oi))**2
def denom(wavelength, TheoSpec,ExperSpec,V_oi):
    return (Y(wavelength, TheoSpec,V_oi)**2) + (Y(wavelength, ExperSpec,V_oi)**2)

def ypendry(TheoSpec,ExperSpec):
    #specDif = SpecDifferenceSq(TheoSpec,ExperSpec)
    return ( ( scipy.integrate.quad(num,0, len(TheoSpec)-1, (TheoSpec,ExperSpec, ImaginaryEnergy(TheoSpec)), limit=100)[0]) / (
    scipy.integrate.quad(denom,0, len(TheoSpec)-1, (TheoSpec,ExperSpec,ImaginaryEnergy(TheoSpec)), limit=100))[0])



def fetchIr(path,column,ran1,ran2):
    logfile =  open(path, 'r')
    logtest = logfile.readlines()

    logfile.close()
    rows = len(logtest[0].split())
    columns = len (logtest)

    IrArray = np.zeros((rows,columns ))
    x= 0
    y = 0
    for line in logtest:


        for word in line.split():
            word = word.replace(',' , '.')
            if (x == 0) and (float(word) < ran1 or float(word) > ran2):

                break


            #print(word, '\n')
            #if float(word) < ran1:
          #      break
           # if float(word) > ran2:
            #    break
            IrArray[x,y] = word
            x+=1
        y+=1
        x=0
    #print(IrArray.shape)
    mask = (IrArray[0,:]) != [0] *len(IrArray[0])
    #print(mask)

    IrArray2 = IrArray[:,mask]
    #print(IrArray2.shape)11
    #plt.scatter(IrArray2[0],IrArray2[1])
    #print(IrArray)

    import matplotlib.image as mpimg
    return IrArray2[(0,column),:]


def IrPlotter(item,title,ran1,ran2,leg = [], multiple = False):
    #colors = np.array([(254,230,206), (253,174,107),(230,85,13)])
    #colors = colors/256
    #colors =['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']
   # print(colors)


    if not(multiple):
        plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),item)
    else:
        index = 0
        for n in item:
            plt.plot(np.linspace(ran1,ran2,(int(ran2-ran1) + 1)),n)
            index += 1
    if len (leg) > 0:
        plt.legend(leg)
    #plt.ylim(bottom=0)
    #plt.gca().invert_yaxis()
    plt.title(title)
    plt.xlabel("cm^-1")
    plt.xlim(ran2,ran1)

    plt.show()
    plt.clf()


#IrPlotter(fetchIr('UntreatedSample.txt',1), "Test")
def gaussian_broadening(spectra, broaden, ran1,ran2,resolution=1,theory=False):

    """ Performs gaussian broadening on IR spectrumimport matplotlib.image as mpimg
    generates attribute self.IR - np.array with dimmension 4000/resolution consisting gaussian-boraden spectrum

    spectra should be in numpy format or list with frequencies in 0 index then intensities in index 1
    :param broaden: (float) gaussian broadening in wn-1"""
    """
    :param resolution: (float) resolution of the spectrum (number of points for 1 wn) defaults is 1, needs to be fixed in plotting
    """
    IR = np.zeros(resolution*(int(ran2-ran1) + 1))
    X = np.linspace(ran1,ran2, resolution*(int(ran2-ran1)+1))

    #for f, i in zip(spectra[:,0]):
      #  IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)
      #  IR=np.vstack((X, IR)).T #tspec

    freq = spectra[0]
    inten = spectra[1]
    if theory:
        freq = spectra[:,0]*.965
        inten = spectra[:,1]



    #print(len(freq))
    for f,i in zip(freq,inten):
       IR += i*np.exp(-0.5*((X-f)/int(broaden))**2)


    return IR
#Define the Gaussian function
def gauss(x, h, A, x0, c):
    return h + (A * np.exp (-1*((x - x0) ** 2 / ( c ** 2))))

IRF = np.zeros((2,(ran2-ran1+1)))
IR1 = fetchIr('UntreatedSample.txt',3,ran1,ran2)
#print('A \n', IR1[1][::10])
IR2 = fetchIr('UntreatedSample.txt',5,ran1,ran2)
IRF [0,:] = gaussian_broadening(IR1,broad,ran1,ran2)

# =============================================================================
# =============================================================================
# for n in range(0,1):
#      IRboy = fetchIr('UntreatedSample.txt',n,ran1,ran2)
#      fit_y = curve_fit(gauss,IRboy[0],IRboy[1],[1,1,1650,1] )
#
#      broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
# #fit_y3 = curve_fit(gaussTwo,IR1[0],IR1[1],p0 = [1,1,1650,1,1,1650,1])
# #print("fit1", fit_y)
# #fit_y =gauss(IR1[0], fit_y[0],fit_y[1],fit_y[2],fit_y[3])\
#      x_range = np.linspace(ran1,ran2,ran2-ran1+1)
#      y_out = gauss (x_range, fit_y[0][0], fit_y[0][1],fit_y[0][2], fit_y[0][3] )
#
# #for x in IR1[0]:
#  #   y_out.append(gauss(x,fit_y[0],fit_y[1],fit_y[2],fit_y[3]))
# #plt.plot(IR1[0],y_out)
#
#      plt.plot(x_range, y_out)#,markersize=1)
#      plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
#      plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, y_out):.5f} ")
#      plt.xlim(max(x_range), min(x_range))
# #plt.xlim(0,4000)
#      plt.title("One Guassian" + " untreated " + str(n))
#      plt.show()
#      plt.clf()
# =============================================================================
     #plt.plot(IRboy[0], IRboy[1])
     #plt.title("untreated" + str(n))
     #plt.show()
     #plt.clf()
# =============================================================================

# =============================================================================
# =============================================================================
# for n in range(0):
#      IRboy = fetchIr('MeOHSample.txt',n,ran1,ran2)
#      fit_y = curve_fit(gauss,IRboy[0],IRboy[1],[1,1,1650,1] )
#      broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
# #fit_y3 = curve_fit(gaussTwo,IR1[0],IR1[1],p0 = [1,1,1650,1,1,1650,1])
# #print("fit1", fit_y)
# #fit_y =gauss(IR1[0], fit_y[0],fit_y[1],fit_y[2],fit_y[3])\
#      x_range = np.linspace(ran1,ran2,ran2-ran1+1)
#      y_out = gauss (x_range, fit_y[0][0], fit_y[0][1],fit_y[0][2], fit_y[0][3] )
#
# #for x in IR1[0]:
#  #   y_out.append(gauss(x,fit_y[0],fit_y[1],fit_y[2],fit_y[3]))
# #plt.plot(IR1[0],y_out)
#
#      plt.plot(x_range, y_out)#,markersize=1)
#      plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
#      plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, y_out):.5f} ")
#      plt.xlim(max(x_range), min(x_range))
# #plt.xlim(0,4000)
#      plt.title("One Guassian" + " MeOH " + str(n))
#      plt.show()
#      plt.clf()
#      #plt.plot(IRboy[0], IRboy[1])
# =============================================================================
     #plt.title("untreated" + str(n))
     #plt.show()
     #plt.clf()
# =============================================================================
#print(np.max(IR1))

for n in range(0):
     IRboy = fetchIr('WA45Sample.txt',n,ran1,ran2)
     fit_y = curve_fit(gauss,IRboy[0],IRboy[1],[1,1,1650,1] )
     broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
#fit_y3 = curve_fit(gaussTwo,IR1[0],IR1[1],p0 = [1,1,1650,1,1,1650,1])
#print("fit1", fit_y)import matplotlib.image as mpimg
#fit_y =gauss(IR1[0], fit_y[0],fit_y[1],fit_y[2],fit_y[3])\
     x_range = np.linspace(ran1,ran2,ran2-ran1+1)
     y_out = gauss (x_range, fit_y[0][0], fit_y[0][1],fit_y[0][2], fit_y[0][3] )

#for x in IR1[0]:
 #   y_out.append(gauss(x,fit_y[0],fit_y[1],fit_y[2],fit_y[3]))
#plt.plot(IR1[0],y_out)

     plt.plot(x_range, y_out)#,markersize=1)import matplotlib.image as mpimg
     plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
     plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, y_out):.5f} ")
     plt.xlim(max(x_range), min(x_range))
#plt.xlim(0,4000)
     plt.title("One Guassian" + " WA45 " + str(n))
     plt.show()
     plt.clf()
     #plt.plot(IRboy[0], IRboy[1])
     #plt.title("untreated" + str(n))
     #plt.show()
     #plt.clf()
# =============================================================================
#print(np.max(IR1))

#=============================================================================
def fourgauss(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2, sigma_2, A_3, x0_3, sigma_3):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))  + \
        ( A_2 * np.exp(-(x - x0_2) ** 2 / ( sigma_2 ** 2)))   + \
    (A_3 * np.exp(-(x - x0_3) ** 2 / (sigma_3 ** 2)))
#=============================================================================
def fivegauss(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1,  A_2, x0_2,  \
              sigma_2, A_3, x0_3, sigma_3, A_4, x0_4, sigma_4 ):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))  + \
        ( A_2 * np.exp(-(x - x0_2) ** 2 / ( sigma_2 ** 2)))   + \
    (A_3 * np.exp(-(x - x0_3) ** 2 / (sigma_3 ** 2))) +   \
       (A_4 * np.exp(-(x - x0_4) ** 2 / (sigma_4 ** 2)))


def gaussTwo(x, H_0, A_0, x0_0, sigma_0,  A_1, x0_1, sigma_1):
    return (H_0 + A_0 * np.exp(-(x - x0_0) ** 2 / (sigma_0 ** 2))) +   \
              ( A_1 * np.exp(-(x - x0_1) ** 2 / (sigma_1 ** 2)))

def fitFive(selec):
    IRboy = fetchIr('MeOHSample.txt',selec,ran1,ran2)
    broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
    print("bob")
    peak1 =  random.randrange(ran1,ran2,1) #1620
    peak2 =  random.randrange(ran1,ran2,1) #1645
    peak3 =   random.randrange(ran1,ran2,1)# # random.randrange(ran1,ran2,1)
    peak4 =  random.randrange(ran1,ran2,1)
    peak5 = 1650 #1678
    print('The peaks are',peak1, peak2,peak3,peak4)
    guesses =[1] +[1,peak1,10]+[1,peak2,10]+[1,peak3,10]+[1,peak4,10]+ [1,peak5,25]
    constraints = ([-20,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5],[20,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,15,np.inf,ran2,50] )

    fit_y5 = curve_fit(fivegauss,IRboy[0],IRboy[1], p0 = guesses,bounds = constraints, method ='trf')
     #print(fit_y2)
    par5 = fit_y5[0]
    yplot5 =fivegauss(x_range, par5[0],par5[1],par5[2],par5[3],par5[4],par5[5],par5[6],par5[7],par5[8],par5[9],par5[10],par5[11], par5[12], par5[13],par5[14], par5[15] )

#print("gauss2",fit_y2)
    plt.plot(x_range,yplot5)
    plt.xlim(max(x_range), min(x_range))
    plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
    plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, yplot5):.5f} ")

     #plt.ylim(bottom=0)
#plt.xlim(0,4000)
#fit_y2 = (curve_fit(fourgauss,IR1[0],IR1[1]))[0]
#plt.plot(IR1[0],fit_y2[1])
#plt.plot()
#plt.xlim(0,4000)
    plt.title('Sum of Five Gaussians')

    plt.show()
    plt.clf()
    gauss1 = gauss(x_range, par5[0], par5[1], par5[2], par5[3])
    gauss2 = gauss(x_range, par5[0], par5[4], par5[5], par5[6])
    gauss3 = gauss(x_range, par5[0], par5[7], par5[8], par5[9])
    gauss4 = gauss(x_range, par5[0], par5[10], par5[11], par5[12])
    gauss5 = gauss(x_range, par5[0], par5[13], par5[14], par5[15])

    gausses = [gauss1,gauss2,gauss3,gauss4,gauss5]
    plt.plot(x_range, gausses[0])
    plt.plot(x_range, gausses[1])
    plt.plot(x_range, gausses[2])
    plt.plot(x_range, gausses[3])
    plt.plot(x_range, gausses[4])
    plt.title("Five Gaussians Methanol")
    peaks0 = find_peaks(gausses[0])[0] +ran1
    peaks1 = find_peaks(gausses[1])[0] +ran1
    peaks2 = find_peaks(gausses[2])[0] +ran1
    peaks3 = find_peaks(gausses[3])[0] +ran1
    peaks4 = find_peaks(gausses[4])[0] +ran1
    plt.legend([f"Gauss 1 {peaks0}cm ",f"Gauss 2 {peaks1}cm ",f"Gauss 3 {peaks2}cm ",f"Gauss 4 {peaks3}cm",f"Gauss 5 {peaks4}cm"] )
    plt.xlabel("cm^-1")
    plt.xlim(max(x_range), min(x_range))
     #plt.ylim(bottom=0)
    plt.show()
    plt.clf()

     #erro = par4[1]
     #print(fit_y2[1])
#for n in range(15):
#fitFive (1)
b = 11
#print(f"{ypendry(IRF [0,:], y_out):.5f}")
Humdities = [5,10,20,30,40,50,60,70,80,90,95]
#for b in range(1):
#    plt.imshow([[1,1,1],[1,1,1],[1,1,1]] )
#    plt.show()
#    plt.clf()
solvents = ["Untreated", "MeOH", "WA45"]
sumsol = np.zeros((4,116))
for sol in range(len(solvents)):
     for hum in range(len(Humdities)):
      IRboy = fetchIr(f'{solvents[sol]}Sample.txt',hum+1,ran1,ran2)
      broadMan = gaussian_broadening(IRboy,broad,ran1,ran2)
      #print("Run", b)
      peak1 =  1620
      peak2 =  1645
      peak3 =   1660# # random.randrange(ran1,ran2,1)
      peak4 =  1678
     # print('The peaks are',peak1, peak2,peak3,peak4)
      guesses =[1] +[1,peak1,10]+[1,peak2,10]+[1,peak3,10]+[1,peak4,10]
      constraints = ([-20,0,ran1,5,0,ran1,5,0,ran1,5,0,ran1,5],[20,np.inf,ran2,15,np.inf,ran2,15,np
                                                                .inf,ran2,15,np.inf,ran2,15] )
      fit_y2 = curve_fit(fourgauss,IRboy[0],IRboy[1], p0 = guesses,bounds = constraints, method ='trf')
      #print(fit_y2)
      par4 = fit_y2[0]
      #erro = par4[1]
      #print(fit_y2[1])
      minSpec = min(fourgauss(x_range, par4[0],par4[1],par4[2],par4[3],par4[4],par4[5],par4[6],par4[7],par4[8],par4[9],par4[10],par4[11], par4[12] ))
      max1 = max (gauss(x_range, par4[0], par4[1], par4[2], par4[3]))
      max2 = max(gauss(x_range, par4[0], par4[4], par4[5], par4[6]))
      max3 = max (gauss(x_range, par4[0], par4[7], par4[8], par4[9]))
      max4 = max(gauss(x_range, par4[0], par4[10], par4[11], par4[12]))
      br = False
      for maxNum in [max1,max2,max3,max4]:
          #print("XD \n", maxNum / (min([max1,max2,max3,max4])-minSpec))
          if maxNum / abs((min([max1,max2,max3,max4])-minSpec)) > 10:
              br = True
      if br:
          print('bad fit')#continue
      if minSpec < 0:
          continue


#print("fit2",fit_y2)
#print("bob", par4)
      yplot4 =fourgauss(x_range, par4[0],par4[1],par4[2],par4[3],par4[4],par4[5],par4[6],par4[7],par4[8],par4[9],par4[10],par4[11], par4[12] )

#print("gauss2",fit_y2)
      plt.plot(x_range,yplot4)
      plt.xlim(max(x_range), min(x_range))
      plt.scatter(IRboy[0],IRboy[1],color='red')#, linewidth = .001 )
      plt.xlabel(f"cm^-1  / Error: {ypendry(broadMan, yplot4):.5f} ")

      #plt.ylim(bottom=0)
#plt.xlim(0,4000)
#fit_y2 = (curve_fit(fourgauss,IR1[0],IR1[1]))[0]
#plt.plot(IR1[0],fit_y2[1])
#plt.plot()
#plt.xlim(0,4000)
      plt.title(f'Sum of Four Gaussians {solvents[sol]} {Humdities[hum]}% humidity')

      plt.show()
      plt.clf()
      #print (max (gauss(x_range, par4[0], par4[1], par4[2], par4[3])))
    # print (max (gauss(x_range, par4[0], par4[4], par4[5], par4[6])))
    # print (max (gauss(x_range, par4[0], par4[7], par4[8], par4[9])))
    # print (max (gauss(x_range, par4[0], par4[10], par4[11], par4[12])))

#OneOf4_gauss =


      peaks4s = ([par4[2], par4[5],par4[8],par4[11]])
      indexL =[]
      for peaky in sorted(peaks4s):
          indexL.append(peaks4s.index(peaky))

      print(peaks4s)

      GaussOrder = np.argsort([par4[2], par4[5],par4[8],par4[11]]    )
      gauss1 = gauss(x_range, par4[0], par4[1], par4[2], par4[3])
      area1 = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[1], par4[2], par4[3] ) )[0]
      #print("The area is " , area1)
      gauss2 = gauss(x_range, par4[0], par4[4], par4[5], par4[6])
      gauss3 = gauss(x_range, par4[0], par4[7], par4[8], par4[9])
      gauss4 = gauss(x_range, par4[0], par4[10], par4[11], par4[12])

      area1 = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[1], par4[2], par4[3] ) )[0]
      area2 = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[4], par4[5], par4[6] ) )[0]
      area3 = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[7], par4[8], par4[9] ) )[0]
      area4 = scipy.integrate.quad(gauss, ran1,ran2, (par4[0], par4[10], par4[11], par4[12] ) )[0]
      print("The areas are " , area1,area2, area3, area4 )
      gausses = [gauss1,gauss2,gauss3,gauss4]
      gausses = [gausses[i] for i in GaussOrder]
      #plt.plot(x_range, gausses[0])
     # plt.plot(x_range, gausses[1])
     # plt.plot(x_range, gausses[2])
      #plt.plot(x_range, gausses[3])
      sumsol [0,:] += gausses[0]
      sumsol [1,:] += gausses[1]
      sumsol [2,:] += gausses[2]
      sumsol [3,:] += gausses[3]




     # plt.title(f"Four Gaussians: ({peaks4s[indexL[0]]:4.1f}, {peaks4s[indexL[1]]:4.1f}, {peaks4s[indexL[2]]:4.1f}, {peaks4s[indexL[3]]:4.1f} {Humdities[b-1]}% humidity)")
      #plt.legend(["Beta Sheet", "Random Coil","Alpha Helix","Beta Turn"])
     # plt.legend([f"Gauss1 {area1:.4}", f"Gauss2 {area2:.4}", f"Gauss3 {area3:.4}", f"Gauss4 {area4:.4}"])
      #plt.xlabel("cm^-1")
      plt.xlim(max(x_range), min(x_range))
      #plt.ylim(bottom=0)
      #plt.show()
     # plt.clf()
      #IRF [0,:] = gaussian_broadening(IR1,broad,ran1,ran2)
      #IRF [1,:] = gaussian_broadening(IR2,broad,ran1,ran2)
      print(f"Error: {ypendry(broadMan, yplot4):.5f}")

     IrPlotter( broadMan, 'Gaussian Broadened', ran1,ran2)
#IrPlotter( IRF [1,:], 'Unstreated Spectra_new_fit', ran1,ran2)
