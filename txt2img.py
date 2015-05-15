# Simply convert all the text files to the image file
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

#dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results47"
dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results46"
os.chdir(dname)
files = os.listdir(dname)



for line in files:
	if (line[-4:] == '.txt' and line[0] == '1'):
		newline = line[:-4]+'.png'
		#if newline not in files:
		namedata = line[0:-4].strip().split('_')
		latt_len = int(namedata[1])
		print latt_len
		lattice = np.loadtxt(line, skiprows = 1,delimiter=',')
		lattice = lattice[0:latt_len,:]
		print lattice.shape
		print newline
		plt.imsave(newline,lattice,[0,2])

plt.show()
