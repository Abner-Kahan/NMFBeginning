#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import random

import re, os
from subprocess import Popen, PIPE
from .utilities import *
import networkx as nx
from operator import itemgetter, attrgetter
import py3Dmol as p3D


class Conformer():

    """
    A class that creates an instance of a molecule defined as conformer.
    It parses gaussian output file (optimization + freq at this moment to
    set proper flags, to be fixed for only freq calcs) and creates an object
    """

    def __init__(self, topol, output_path):
        """Construct a conformer object
        :param topol: 
        :param output_path: (string) this specifies the directory any generated IR plots will be placed
        """
        self._id = topol
        self.topol = topol
        self.output_path = output_path

    def load_model(self, file_path):
        """Loads a file and constructs a conformer
        :param file_path: (string) this is the directory of the data file being used to condtruct a conformer object
        """
        self.NAtoms = None
        self._id    = str(file_path).split('/')[-2]
        self.topol = self._id
        geom = [] ; atoms = []

        for n, line in enumerate(open(file_path, 'r').readlines()):

            if n == 0 and self.NAtoms == None: self.NAtoms = int(line)
            if n > 1:
                if len(line.split()) == 0: break 
                geom.append([float(x) for x in line.split()[1:4]])
                if line.split()[0].isalpha(): atoms.append(line.split()[0])
                else: atoms.append(element_symbol(line.split()[0]))

        self.xyz = np.array(geom)
        self.atoms = atoms
    def load_log(self, file_path):
            """ Creates a conformer object using infromation from a file given the file path
            :param file_path: (string) path to file
            """
            # why did this try function get commented out?
            #try:
            #    logfile = open(file_path, 'r')
            #except IOError: 
            #    print("%30s not accessible", file_path)
            #    return 1 
                
            normal_mode_flag=False
            freq_flag = False
            read_geom = False
            opt_flag = False
    
            #temprorary variables to hold the data
            freq = [] ; ints = [] ; vibs = [] ; geom = [] ; atoms = []
    
            job_opt = False ; job_freq = False ; job_optfreq = False ; job_sp = False ; job = 0
    
            self.NAtoms = None
            self._id    = file_path.split('/')[-2]
            self.path   = '/'.join(file_path.split('/')[:-1])
    
            for line in open(file_path, 'r').readlines():
    
                    if re.search('^ #', line) and job == 0:
                        if "opt" in line:
                            if "freq" in line: 
                                job_optfreq = True ; job += 1 
                                #print("Reading optfreq")
                            else: 
                                job_opt = True ; job += 1 
                                #print("Reading opt")
                        elif "freq" in line: 
                                job_optfreq = True ; freq_flag = True ;  job += 1 
                                job_freq = True
                                #print("Reading freq")
                        else: job_sp = True ; job += 1 
    
                    if self.NAtoms is None and re.search('^ NAtoms=', line): 
                        self.NAtoms = int(line.split()[1])
                        self.NVibs  = self.NAtoms*3-6
    
                    if self.NAtoms is None and job_freq == True and re.search('Deg. of freedom', line):
                        self.NVibs  = int(line.split()[3])
                        self.NAtoms = int((self.NVibs + 6)/3)
    
                    if re.search('^ Frequencies', line):        
                        freq_line = line.strip() 
                        for f in freq_line.split()[2:5]: freq.append(float(f))
                        normal_mode_flag = False
    
                    elif re.search('^ IR Inten', line):        
                        ir_line = line.strip()                 
                        for i in ir_line.split()[3:6]: ints.append(float(i))
    
                    elif re.search('^  Atom  AN', line): 
                         normal_mode_flag = True          #locating normal modes of a frequency
                         mode_1 = []; mode_2 = []; mode_3 = [];
                         continue
    
                    elif normal_mode_flag == True and re.search('^\s*\d*\s*.\d*', line) and len(line.split()) > 3:
    
                         mode_1.append([float(x) for x in line.split()[2:5]])
                         mode_2.append([float(x) for x in line.split()[5:8]])
                         mode_3.append([float(x) for x in line.split()[8:11]])
    
                    elif normal_mode_flag == True: 
    
                         normal_mode_flag = False 
                         for m in [mode_1, mode_2, mode_3]: vibs.append(np.array(m))
    
    
                    if job_optfreq == True:
    
                        if freq_flag == False and re.search('Normal termination', line): freq_flag = True
    
                        elif freq_flag == True and re.search('SCF Done',   line): self.E = float(line.split()[4]) 
                        elif freq_flag == True and re.search('Sum of electronic and zero-point Energies',   line): self.Ezpe = float(line.split()[6])
                        elif freq_flag == True and re.search('Sum of electronic and thermal Enthalpies' ,   line): self.H    = float(line.split()[6])                    
                        elif freq_flag == True and re.search('Sum of electronic and thermal Free Energies', line): self.F    = float(line.split()[7])
    
                        elif freq_flag == True and re.search('Coordinates', line): 
                            if len(geom) == 0: read_geom = True
                        elif freq_flag == True and read_geom == True and re.search('^\s*.\d', line):
                            #geom.append(map(float, line.split()[3:6])) 
                            #convert to a parse directly into list rather than map
                            geom.append([float(x) for x in line.split()[3:6]])
                            atoms.append(element_symbol(line.split()[1]))
                            if int(line.split()[0]) == self.NAtoms:
                               read_geom = False
    
                    elif job_opt == True: 
    
                        if re.search('SCF Done',   line): E = float(line.split()[4])
                        if re.search('Optimization completed.', line): 
                             self.E = E ; opt_flag = True  
                        elif opt_flag == True and re.search('Coordinates', line) : read_geom = True
                        elif opt_flag == True and read_geom == True and re.search('^\s*.\d', line):
                             #geom.append(map(float, line.split()[3:6])) 
                             #convert to a parse directly into list rather than map
                             geom.append([float(x) for x in line.split()[3:6]])
                             atoms.append(element_symbol(line.split()[1]))
                             if int(line.split()[0]) == self.NAtoms:
                               read_geom = False
    
                    elif job_sp == True:
    
                        print("No idea what you're dong")
    
            if freq_flag == True: 
                self.Freq = np.array( freq ) 
                self.Ints = np.array( ints )
                self.Vibs=np.zeros((self.NVibs, self.NAtoms, 3))
                for i in range(self.NVibs): self.Vibs[i,:,:] = vibs[i]
    
            self.xyz = np.array(geom)
            self.atoms = atoms
     
Conformer.load_log(Conformer,'Tri_A131/Tri_A131_0000/input.log')