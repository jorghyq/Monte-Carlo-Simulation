# Program to run the
import os
import sys
import subprocess
import shlex
#from subprocess import call
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time

total_run = 1000000
num_mol = 100
num_metal = 0
cenergy = 20
venergy = 5
mcenergy = 20
ffn=100
latt_len = 100

if len(sys.argv) > 1:
    total_run = float(sys.argv[1])
elif len(sys.argv > 2):
    total_run = float(sys.argv[1])
    num_metal = float(sys.argv[2])
command = '/home/jorghyq/Project/Monte-Carlo-Simulation/mc-rect-lattice-func-linux4 -a %d -b %d -c %d -d %d -e %f -f %d -g %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy,ffn)
print command
args = shlex.split(command)
print args
#p = subprocess.check_call(args,shell=True)
os.system(command)
path = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results100/"
fname = os.listdir(path)
print fname,fname[-1]
lattice = np.loadtxt(path+fname[-1],delimiter=',',skiprows=1)
#lattice = np.zeros((100,100))
#lattice = np.loadtxt("D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results19\\1.0e+009_100_400_250_4.0e+001_5.0e+000_4.0e+001.txt", delimiter=',',skiprows=1)
#lattice = lattice[0:latt_len,0:]
temp1, temp2 = np.where(lattice == 6) # mol1
temp3, temp4 = np.where(lattice == 4) # mol2
temp5, temp6 = np.where(lattice == 1)  # metal
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

ax.scatter(temp1,temp2,s = 50,c = "#0ACEF5",linewidth='0',marker = xy1)
ax.scatter(temp3,temp4,s = 50,c = "b",linewidth='0',marker = xy1)
ax.scatter(temp5,temp6,s = 5,c = "#F78C00",linewidth='0',marker = "o")
plt.xlim([0, 100])
plt.ylim([0, 100])

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
#ax.text(-10,0,headdate[0])
ax.set_aspect("equal")
plt.show()
