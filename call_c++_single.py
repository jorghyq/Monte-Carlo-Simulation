# Program to run the
import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time

total_run = 10000000
latt_len = 80
num_mol = 30
num_metal = 30
cenergy = 24
venergy = 0
mcenergy = 24
os.system('mc-rect-lattice-func -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
lattice = np.loadtxt("latt.txt", delimiter=',')
temp1, temp2 = np.where(lattice == 2)
temp3, temp4 = np.where(lattice == 1)
plt.figure()
plt.imshow(lattice)
#print line
#fig = plt.figure()
#ax = plt.axes()
#fig.add_axes(ax)
#ax.set_aspect("equal")
#ax.scatter(temp1,temp2,s = 70,c = "r",marker = "s")
#ax.scatter(temp3,temp4,s = 60,c = "b",marker = "s")
#plt.show()
		#newline = line[:-4]+'.png'
		#plt.savefig(newline)
		#print newline
		#plt.imsave(newline,lattice,[0,2])
plt.show()
