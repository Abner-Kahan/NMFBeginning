BGH2-difference plots the scaled and actual difference between the high and low portions of the IR spectra for the BGH2 silk model.

BGH2NMF plots the scaled difference between the high and low portions of the IR spectra for the BGH2 silk model, as well as decompose the higher frequency portion of the spectrum into two components using NMF.
The user can adjust the parameters of the NMF.

Gaussian8Single.py first plots an Gaussian broadened spectra on the user's choice in humidity (0:5, 1:10, 2:20, 3:30, 4:40, 5:50, 6:60, 7:70, 8:80, 9:90, 10:95) and the solvent (0:untreated, 1:methanol, 2:water), which the user can control. This plot is over the range of both Amide I and Amide II regions.
It then plots an overlaid graph including red dots representing a fit using the sum of 8 Gaussians.
It then plots the individual 8 Gaussians in the next figure, as well as the sums of the first four Gaussians for the Amide I band and the next four Gaussians for the Amide II band. The program also returns the wavenumber of all the found peaks.

Gaussian4AmideISingle.py returns red dot representing the sum of four fitted Gaussians for one sample overlaid on top of a Gaussian broadened Amide I region. The user can choose which untreated spectra to analyze or change the file name get wa45 and methanol treated. The ypendry error is also on the figure. The next figure outputted plots each of the individual 4 gaussians, their wavenumber, and their area.

Gaussian4AmideIAll 
