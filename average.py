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


totrun()
#def aver():



