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
import json

def sorted_ls(path, files):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(files, key=mtime))

total_run = 1000000
num_mol1 = 100
num_mol2 = 0
num_metal = 0
cenergy = 40
venergy = 5
mcenergy = 40
ffn=100
latt_len = 100
restore = 1
print len(sys.argv)
for i in range(len(sys.argv)):
    print i, sys.argv[i]
if len(sys.argv) == 2:
    total_run = float(sys.argv[1])
elif len(sys.argv) ==  3:
    total_run = float(sys.argv[1])
    num_mol1 = float(sys.argv[2])
elif len(sys.argv) == 4:
    total_run = float(sys.argv[1])
    num_mol1 = float(sys.argv[2])
    num_mol2 = float(sys.argv[3])
elif len(sys.argv) == 5:
    total_run = float(sys.argv[1])
    num_mol1 = float(sys.argv[2])
    num_mol2 = float(sys.argv[3])
    num_metal = float(sys.argv[4])
elif len(sys.argv) == 6:
    total_run = float(sys.argv[1])
    num_mol1 = float(sys.argv[2])
    num_mol2 = float(sys.argv[3])
    num_metal = float(sys.argv[4])
    ffn = int(sys.argv[5])
else:
    pass

print total_run, num_mol1, num_mol2, num_metal
command = './mc-rect-lattice-func-linux6 -a %d -b %d -c %d -d %d -e %f -f %f -g %f -h %d -i %d' % (total_run,num_mol1,num_mol2,num_metal,cenergy,venergy,mcenergy,ffn,restore)
args = shlex.split(command)
print args
#p = subprocess.check_call(args,shell=True)
os.system(command)
fpath = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results%d/" % ffn
files = os.listdir(fpath)
print fpath
fname = sorted_ls(fpath,files)
if fname[-2][-5] == 't':
    print fname[-1]
    lattice = np.loadtxt(fpath+fname[-1],delimiter=',',skiprows=1)
else:
    print fname[-2]
    lattice = np.loadtxt(fpath+fname[-2],delimiter=',',skiprows=1)
#lattice = np.zeros((100,100))
#lattice = np.loadtxt("D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results19\\1.0e+009_100_400_250_4.0e+001_5.0e+000_4.0e+001.txt", delimiter=',',skiprows=1)
# To have a customer designed shape, one has to draw it by himself
x = [-1.5,-0.5,-0.5,0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5]
y = [0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5,-0.5,-0.5,0.5]
x2 = [-3,-1,-1,1,1,3,3,1,1,-1,-1,-3,-3]
y2 = [1,1,3,3,1,1,-1,-1,-3,-3,-1,-1,1]
x3 = [-0.5,0.5,0.5,-0.5,-0.5]
y3 = [0.5,0.5,-0.5,-0.5,-0.5]
#x3 = x/2
#y3 = y/2
xy1 = list(zip(x,y))
xy2 = list(zip(x2,y2))
xy3 = list(zip(x3,y3))

# load colors
with open('colors.json','r') as f:
    colors = json.load(f)

fig = plt.figure()
ax = plt.axes()
fig.add_axes(ax)


# plot metal
coor_x, coor_y = np.where(lattice == 1) # metal
ax.scatter(coor_x,coor_y,s=20,c=colors['1'],linewidth='0',marker = "o")

# plot molecular center
for i in [11,12]:
    coor_x,coor_y = np.where(lattice == i)
    ax.scatter(coor_x,coor_y,s=120,c=colors[str(i)],linewidth='0',marker=xy1)

# plot endgroups
for i in range(2,11):
    coor_x,coor_y = np.where(lattice == i)
    if str(i) in colors:
        ax.scatter(coor_x,coor_y,s=20,c=colors[str(i)],linewidth='0',marker=xy3)

plt.xlim([-0.5, latt_len-0.5])
plt.ylim([-0.5, latt_len-0.5])

ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
#ax.text(-10,0,headdate[0])
ax.set_aspect("equal")
plt.show()
