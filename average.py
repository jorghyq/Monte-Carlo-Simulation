# average
import os
from analyzer import McAnalyzer
import platform
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster
from PIL import Image

filenumbers = [15,16,18,27,28,29,30,31,32,33,34,35,36,42,43,44,45]

temp_path = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results"

#totalenergy = np.zeros((9,13))
#totalenergy_av = np.zeros((9,13))
#mdense = np.zeros((9,13))
#m1d = np.zeros((9,13))
#m2d = np.zeros((9,13))
#mdis = np.zeros((9,13


# calculate the results
def totrun():
    for temp_file in filenumbers:
        az = McAnalyzer('None')
        full_path = temp_path + str(temp_file)
        print full_path
        az.set_path(full_path)
        print az.path
        az.load_logfile()
        az.run(1)


#totrun()
class McAverager:
    def __init__(self,wd_path,analyzer):
        self.path = wd_path
        self.analyzer = analyzer
        #self.analyzer.set_path(self.path)
        #self.analyzer.load_logifle()

        # read all the files
        #files = os.listdir(self.path)
        # select the files
        #for file in files:
        #    if file == "logfile.txt":



    def set_path(self, new_path):
        self.path = new_path

    def set_filenumbers(self,filenumbers):
        self.filenumbers = filenumbers

    def init_range(self):
        os.chdir(self.path)
        self.mdense = np.loadtxt("mdense.txt", delimiter=',')
        self.m1d = np.loadtxt("m1d.txt", delimiter=',')
        self.m2d = np.loadtxt("m2d.txt", delimiter=',')
        self.mdis = np.loadtxt("mdis.txt", delimiter=',')
        self.totalenergy = np.loadtxt("totalenergy.txt",delimiter=',')
        self.totalenergy_av = np.loadtxt("totalenergy_av.txt",delimiter=',')
        self.metal_ind = (self.mdense.shape)[0]
        self.energy_ind = (self.mdense.shape)[1]

    def averg(self):
        for filenum in filenumbers:
            fulldir = self.path + str(filenum)
            os.chdir(fulldir)
            temp_mdense = np.loadtxt("mdense.txt", delimiter=',')
            temp_m1d = np.loadtxt("m1d.txt", delimiter=',')
            temp_m2d = np.loadtxt("m2d.txt", delimiter=',')
            temp_mdis = np.loadtxt("mdis.txt", delimiter=',')
            temp_totalenergy = np.loadtxt("totalenergy.txt",delimiter=',')
            temp_totalenergy_av = np.loadtxt("totalenergy_av.txt",delimiter=',')
            self.mdense = self.mdense + temp_








