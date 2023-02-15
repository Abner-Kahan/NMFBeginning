SingleTesterBinaryMultiple generates .npy files that contain decompositions for various sugars into various components. One can adjust the sugarList to apply NMF decomposition to a different range of sugars, each of which needs to have their number of residues specified in residueCount List.

[A full list of available sugars are found in the VMD subdirectory which contains distance maps for sugars with 1's being represented if two residues are within 4.5 Angstrom of another residue. Those are labeled as 45d.txt]
One could also adjust the number of components tested for the decompositions, in the list for which k iterates through. NMF hyper parameters can be adjusted in the NMFmap function.
This program generates labeled .npy files containing the word ProCompons to represent the residues maps of the glycans, and TempCompons to represent the degree to which the interactions fluctuate over time. The name of the spy files also contains the sugar and the number of NMF components.
Examples of the decompositions are found in MD decomposs subdirectory in this folder.

AnalyzeMultiple