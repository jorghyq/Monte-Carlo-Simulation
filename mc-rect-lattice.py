# Monte Carlo

import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import copy
import time
#import Image

#global lattice_size
LATTICE_CONSTANT = 1
TOTAL_RUN = 10000
print "Total run: %d" % TOTAL_RUN

def set_in_range(x,range):
    if x > range-1:
        x = x - range
    elif x < 0:
        x = x + range
    return x


def get_coor_mol(x,y):
    coor = np.zeros((5,2))
    coor[0,:] = [x, y]
    coor[1,:] = [set_in_range(x-1,latt_len), y]
    coor[2,:] = [set_in_range(x+1,latt_len), y]
    coor[3,:] = [x, set_in_range(y-1,latt_len)]
    coor[4,:] = [x, set_in_range(y+1,latt_len)]
    return coor

def set_mol(x,y,op,id,coor,latt):
    temp = get_coor_mol(x,y)
    for i in range(0,len(temp)):
        latt[temp[i][0],temp[i][1]] = op
    coor[id,:] = temp[0,:]


def is_occupied(x,y,latt):
    temp = get_coor_mol(x,y)
    for i in range(0,len(temp)):
        if latt[temp[i][0],temp[i][1]] == 1:
            return True
    return False
################################### define the lattice ###################
# assume the coordinate is (i,j)
# for a square lattice, points around it is : (i+1,j), (i,j+1), (i-1,j), (i,j-1)
# for a hexagonal lattice, points around it is : (i-1,j+1), (i,j+1), (i+1,j), (i+1,j-1),(i,j-1),(i-1,j)
#
####### How to achieve periodic condition
# for a square lattice: if i < 0, i = i + latticeSize, if i > latticeSize-1, i = i -latticeSize
# for a hexagonal lattice: if i < 0, i = i + latticeSize, if i > latticeSize-1, i = i -latticeSize
#
##########################################################################
latt_len = 40
lattice = np.zeros((latt_len,latt_len))
num_mol = 40
num_metal = 0
coor_mol = np.zeros((num_mol,2))
################### Distribute the molecules ####################
for i in range(0,num_mol):
    state = True
    while state == True:
        ind_x = rd.randint(0, latt_len-1)
        ind_y = rd.randint(0, latt_len-1)
        if is_occupied(ind_x,ind_y,lattice) == False:
            set_mol(ind_x,ind_y,1,i,coor_mol,lattice)
            state = False

print "molecules are distributed"
count = 0
while count < TOTAL_RUN:
    ind_element = rd.randint(0, num_mol)
    energy_current = cal_energy_mol(coor_mol[ind_element,:])






plt.figure()
plt.imshow(lattice)
plt.show()





