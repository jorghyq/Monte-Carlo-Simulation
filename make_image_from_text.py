# Make images from the text files
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results2"
os.chdir(dname)
#print dname
#files  = os.listdir(dname)
#print files
for line in os.listdir(dname):
	if (line[-4:] == '.txt'):
		lattice = np.loadtxt(line, delimiter=',')
		temp1, temp2 = np.where(lattice == 2)
		temp3, temp4 = np.where(lattice == 1)
		#print line
		fig = plt.figure()
		ax = plt.axes()
		fig.add_axes(ax)
		ax.set_aspect("equal")
		ax.scatter(temp1,temp2,s = 70,c = "r",marker = "s")
		ax.scatter(temp3,temp4,s = 60,c = "b",marker = "s")
		newline = line[:-4]+'.png'
		plt.savefig(newline)
		print newline
		plt.imsave(newline,lattice,[0,2])
	
