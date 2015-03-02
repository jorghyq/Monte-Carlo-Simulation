# average
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster
from PIL import Image

filenumbers = [15,16,18,27,28,29,30,31,33]

temp_path = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results"

totalenergy = np.zeros((9,13))
totalenergy_av = np.zeros((9,13))
mdense = np.zeros((9,13))
m1d = np.zeros((9,13))
m2d = np.zeros((9,13))
mdis = np.zeros((9,13

for temp_file in filenumbers:
	full_path = temp_path + str(temp_file)
	
