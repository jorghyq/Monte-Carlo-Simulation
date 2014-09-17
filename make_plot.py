# make plot from data
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results6"
os.chdir(dname)

files = os.listdir(dname)
xvalue = []
yvalue = []
for line in files:
	if (line[-4:] == '.txt'):
		
		namedata = line[0:-4].strip().split('-')
		cenergy = namedata[4]
		print namedata
		print cenergy
		f = open(line, 'r')
		headdata = f.readline(1).strip().split(',')
		latt_len = int(headdata[-1])
		print headdata
		f.close()
