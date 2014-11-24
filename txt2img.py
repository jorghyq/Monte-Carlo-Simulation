# Simply convert all the text files to the image file
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results16\\4-fold-ev-3"
os.chdir(dname)
files = os.listdir(dname)



for line in files:
	if (line[-4:] == '.txt'):
		newline = line[:-4]+'.png'
		#if newline not in files:
		namedata = line[0:-4].strip().split('-')
		latt_len = int(namedata[1])
		print latt_len
		lattice = np.loadtxt(line, skiprows = 1,delimiter=',')
		lattice = lattice[0:latt_len,:]
		print lattice.shape
		print newline
		plt.imsave(newline,lattice,[0,2])

plt.show()
