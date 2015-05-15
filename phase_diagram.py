# make the phase transition diagram

import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
from PIL import Image

#dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results10"
dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results15"
os.chdir(dname)
mdense = np.loadtxt("mdense.txt", delimiter=',')
m1d = np.loadtxt("m1d.txt", delimiter=',')
m2d = np.loadtxt("m2d.txt", delimiter=',')
mdis = np.loadtxt("mdis.txt", delimiter=',')


nmetal_ind = mdense.shape[0]
cenergy_ind = mdense.shape[1]
print nmetal_ind
print cenergy_ind

mtotal = mdense + m1d + m2d + mdis
#print mtotal
# left and right


mdense_p = mdense/mtotal
m1d_p = m1d/mtotal
m2d_p = m2d/mtotal
# up down
mdense_p = np.flipud(mdense/mtotal)
m1d_p = np.flipud(m1d/mtotal)
m2d_p = np.flipud(m2d/mtotal)
#mdense_p = mdense_p[:,:-1]
#m1d_p = m1d_p[:,:-1]
#m2d_p = m2d_p[:,:-1]

#mdense_p = np.fliplr(mdense_p)
#m1d_p = np.fliplr(m1d_p)
#m2d_p = np.fliplr(m2d_p)
#mdense_p = np.transpose(mdense/mtotal)
#m1d_p = np.transpose(m1d/mtotal)
#m2d_p = np.transpose(m2d/mtotal)
r = m1d_p
g = m2d_p
b = mdense_p
#mdis_p = mdis/mtotal
#colors = np.zeros((nmetal_ind,cenergy_ind,3))
rgb = np.dstack((r,g,b))
#print rgb
#print rgb.shape
#print rgb
im = Image.fromarray(np.uint8(rgb*255.999))
#print im
plt.figure()
#plt.scatter()
plt.imshow(im,extent=[0.125,1,0,2],aspect="auto")
plt.xlabel("Eo/Ec")
plt.ylabel("Metal/Molecule")
img = "BDS Phase full.png"
plt.savefig(img,dpi=1000)
plt.show()
