# Program to run the
import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time

total_run = 100000000

num_mol = 400
num_metal = 500
cenergy = 40
venergy = 5
mcenergy = 40
latt_len = 100
#f = open('D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results14\\1.0e+007-80-200-250-13-3.0e+000-13.txt', 'r')
#temp = f.readline().strip()
#headdate = temp.split(',')
#latt_len = int(headdate[-1])
#print headdate
#f.close()
os.system('mc-rect-lattice-func3 -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
lattice = np.loadtxt("latt.txt", delimiter=',',skiprows=1)
lattice1 = lattice[0:latt_len,0:]
lattice2 = lattice[latt_len:,0:]
temp1, temp2 = np.where(lattice1 == 5)
temp3, temp4 = np.where(lattice1 == 1)
temp5, temp6 = np.where(lattice1 == 3)
plt.imshow(lattice1)
plt.imshow(lattice2)
# To have a customer designed shape, one has to draw it by himself
x = [-1.5,-0.5,-0.5,0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5]
y = [0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5,-0.5,-0.5,0.5]
xy1 = list(zip(x,y))
fig = plt.figure()
ax = plt.axes()
fig.add_axes(ax)

ax.scatter(temp1,temp2,s = int(2000/latt_len),c = "r",marker = (xy1,0))
ax.scatter(temp5,temp6,s = int(2000/latt_len),c = "y",marker = (xy1,0))
ax.scatter(temp3,temp4,s = int(1000/latt_len),c = "b",marker = "s")
#ax.text(-10,0,headdate[0])
ax.set_aspect("equal")
plt.show()
