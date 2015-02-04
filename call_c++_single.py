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

num_mol = 200
num_metal = 200
cenergy = 20
venergy = 5
mcenergy = 20
latt_len = 100
#os.system('mc-rect-lattice-func -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
lattice = np.loadtxt("D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results19\\1.0e+009_100_400_250_4.0e+001_5.0e+000_4.0e+001.txt", delimiter=',',skiprows=1)
lattice = lattice[0:latt_len,0:]
temp1, temp2 = np.where(lattice == 3)
temp3, temp4 = np.where(lattice == 1)
# To have a customer designed shape, one has to draw it by himself
x = [-1.5,-0.5,-0.5,0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5]
y = [0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5,-0.5,-0.5,0.5]
x2 = [-3,-1,-1,1,1,3,3,1,1,-1,-1,-3,-3]
y2 = [1,1,3,3,1,1,-1,-1,-3,-3,-1,-1,1]
#x3 = x/2
#y3 = y/2
xy1 = list(zip(x,y))
xy2 = list(zip(x2,y2))
fig = plt.figure()
ax = plt.axes()
fig.add_axes(ax)

ax.scatter(temp1,temp2,s = 20,c = "#0ACEF5",linewidth='0',marker = xy2)
#ax.scatter(temp1,temp2,s = 20,c = "r",linewidth='0',marker = xy2)
ax.scatter(temp3,temp4,s = 5,c = "#F78C00",linewidth='0',marker = "o")
plt.xlim([0, 100])
plt.ylim([0, 100])

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
#ax.text(-10,0,headdate[0])
ax.set_aspect("equal")
plt.show()
