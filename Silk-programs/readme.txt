BGH2-difference plots the scaled and actual difference between the high and low portions of the IR spectra for the BGH2 silk model.

BGH2NMF plots the scaled difference between the high and low portions of the IR spectra for the BGH2 silk model, as well as decompose the higher frequency portion of the spectrum into two components using NMF.
The user can adjust the parameters of the NMF.

Gaussian8Single.py first plots an Gaussian broadened spectra on the user's choice in humidity (0:5, 1:10, 2:20, 3:30, 4:40, 5:50, 6:60, 7:70, 8:80, 9:90, 10:95) and the solvent (0:untreated, 1:methanol, 2:water), which the user can control. This plot is over the range of both Amide I and Amide II regions. Users can adjust the window.
It then plots a graph including red dots representing the IR data points, an an overlaid blue line representing the sum of 8 Gaussians.
It then plots the individual 8 Gaussians in the next figure, as well as the sums of the first four Gaussians for the Amide I band and the next four Gaussians for the Amide II band. The program also returns the wavenumber of all the found peaks.

Gaussian4AmideISingle.py returns a blue line representing the sum of four fitted Gaussians for one sample overlaid on red dots representing the red dated. The user can choose which humidity of untreated spectra to analyze or change the file name to wa45 or methanol treated to change the solvent. The ypendry error is also on the figure. The next figure outputted plots each of the individual 4 gaussians, their wavenumber, and their area. Users can adjust the window.

Gaussian4AmideIAll iterates through all combinations of solvation and humidities started with untreated, the the aforementioned humidities of 5-95, and repeats the process for all the solvent conditions. For each combination, a figure is generated for the Amide I region showing a blue line representing the representing the sum of four Gaussians, an an overlaid set of red dots representing the actual data. On the bottom of the figure the ypendry error is displayed. Users can adjust the window, as well as the initialization values for the peaks. The area of and wavenumber of the four Gaussians is printed in the terminal.

Gaussian9All iterates through all combinations of solvation and humidities started with untreated, the the aforementioned humidities of 5-95, and repeats the process for all the solvent conditions.
For each combination, a figure is generated for the combined Amide I and Amide II region showing a blue line representing the sum of 9 Gaussians, and underlaid red dots representing the IR data points. The ypendry error is then displayed at the bottom of the figure.
After each figure of the sum, an additional figure is shown graphing the nine individual 9 Gaussians. The first four are added together to make the Amide I band, which is displayed on top of those four. This fifth Gaussian is labeled as a correction. The Amide II band is also displayed which is the sum of the four remaining spectra. The wavelengths and areas of all 9 peaks are printed from smallest to largest wavenumber. At the end of the program, a table is printed with the very peak location for the nine gaussians for each solvent conditions, which is saved (with an erroneous name) in Peaksof8A.csv, while the areas for the nine gaussians are  printed and saved into Areasof8A.csv

Silk5NMF plots all the spectra from 0 to 4000 cm-1 across all 33 combinations of solvent conditions.
It then performs NMF decomposition across this window into 5 constituent spectra across the range and displays each of the five constituents individually. One could adjust NMF hyper parameters, the window, and even oscillate through a smaller number of spectra. To adjust the number of components, one adjusts the line when making the NMF model, and displays more or less peak plots.


NMF_AllSilk.py displays two figures. The first figure shows the Amide I region of the Silk for all samples across every condition of humidity. It then shows an NMF decomposition of the amide I band into four labeled components along with an overlaid example of a spectra being decomposed. The wavenumber of the peaks are printed. The hyper parameters of NMF and window can be adjusted.
The next figure shows the percentage of each secondary structure across every solvent and humidity condition. One could correlate the peaks in the last figure, to the percentages in this figure.

NMF_SomeSilk.py  performs the same functions as allSilkApril.py except it easier to select the spectra that are used for the decomposition.
At the beginning of nmf2TesterMixB one changes the first line representing the total number of input spectra.
One then comments and then uncomments the lines representing which spectra are being imported, whether untreated, methanol, or WA45. One needs to make sure that the index numbers to IRF, the import number are correct. One then makes sure the number of subplots is correct and corresponds to one more then the number of solvents being tested (4 for untreated, MeOH, and WA45). One the sets the subplots to plot the corresponding entries of IRF.
Finally one needs to change the fig2, axs2 to get the correct number of subplots corresponding to the solvent conditions being tested, and then graphs the corresponding entries of H


Simulated silk spectra
