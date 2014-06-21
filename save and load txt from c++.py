# python program to load the file and display the file from c++
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import copy
import time

lattice_size = 50
LATTICE_CONSTANT = 1
lattx = np.zeros((lattice_size, lattice_size))
latty = np.zeros((lattice_size, lattice_size))
tempx = 0
tempy = 0
for i in range(0, lattice_size):
	tempx = LATTICE_CONSTANT * math.cos(3.14/3) * i
	for j in range(0, lattice_size):
		lattx[i, j] = tempx
		latty[i, j] = tempy
		tempx = tempx + LATTICE_CONSTANT
	tempy = tempy - LATTICE_CONSTANT * math.sin(3.14/3)	


lattice = np.loadtxt("lat.txt", delimiter=',')
plt.figure()
for i in range(0, lattice_size):
	plt.scatter(lattx[i,:], latty[i,:])

temp1, temp2 = np.where(lattice == 3)
temp3, temp4 = np.where(lattice == 2)
temp5, temp6 = np.where(lattice == 1)
print len(temp1)
plt.scatter(lattx[temp1,temp2],latty[temp1,temp2],s = 70,c = "r",marker = "s")
plt.scatter(lattx[temp3,temp4],latty[temp3,temp4],s = 60,c = "y",marker = "s")
plt.scatter(lattx[temp5,temp6],latty[temp5,temp6],s = 60,c = "k",marker = "s")
#im = plt.imshow(lattice)
plt.show()
