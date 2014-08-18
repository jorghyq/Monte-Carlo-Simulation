# Monte Carlo

import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import copy
import time
# import Image

# global lattice_size
LATTICE_CONSTANT = 1
TOTAL_RUN = 1000000
print "Total run: %d" % TOTAL_RUN

def set_in_range(x,range):
    if x > range-1:
        x = x - range
    elif x < 0:
        x = x + range
    return x


def get_coor_mol(input_coor):
    x = input_coor[0]
    y = input_coor[1]
    coor = np.zeros((5,2))
    coor[0,:] = [x, y]
    coor[1,:] = [set_in_range(x-1,latt_len), y]
    coor[2,:] = [set_in_range(x+1,latt_len), y]
    coor[3,:] = [x, set_in_range(y-1,latt_len)]
    coor[4,:] = [x, set_in_range(y+1,latt_len)]
    return coor

def set_mol(input_coor,op,id_mol,coor,latt):
    temp = get_coor_mol(input_coor)
    for i in range(0,len(temp)):
        latt[temp[i][0],temp[i][1]] = op
    if input_coor[0] != -1:
        coor[id_mol,:] = temp[0,:]

def set_element(input_coor,op,id_ele,coor,latt):
    ele_length = len(input_coor)
    for i in range(0,ele_length):
        latt[input_coor[i][0],input_coor[i][1]] = op
    coor[id_ele,:] = input_coor[0,:]


def set_metal(input_coor,op,id_metal,coor,latt):
    temp = input_coor
    latt[temp[0],temp[1]] = op
    coor[id_metal,:] = temp


def is_occupied(input_coor,latt):
    ele_length = len(input_coor)
    if ele_length == 1:
        latt[input_coor]
    for i in range(0,ele_length):
        if latt[input_coor[i][0],input_coor[i][1]] != 0:
            return True
    return False


def cal_energy_mol(coor,coor_mol,latt):
    # First detect if there are molecules around
    energy = 0
    temp = get_coor_mol(coor)
    # Get the points around this molecule
    pos_around = np.zeros((8,2))
    pos_around[0:3,:] = get_coor_mol(temp[3,:])[1:4,:]
    pos_around[3:5,:] = get_coor_mol(temp[1,:])[[1,4],:]
    pos_around[5:7,:] = get_coor_mol(temp[2,:])[[2,4],:]
    pos_around[7,:] = get_coor_mol(temp[4,:])[4,:]
    mol_around = []
    for i in range(0,8):
        if lattice_num[pos_around[i,0],pos_around[i,1]] != 0:
            mol_around.append(lattice_num[pos_around[i,0],pos_around[i,1]])
    mol_around = list(set(mol_around))
    for i in mol_around:
        temp_coor = coor_mol[i,:]
        dis = math.sqrt(math.pow(temp_coor[0] - coor[0],2) +math.pow(temp_coor[1] - coor[1],2))
        if dis < 2.5:
            energy = energy - 5
        elif dis == 3:
            energy = energy - 2.5
    return energy

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
lattice_num = np.zeros((latt_len,latt_len))

num_mol = 40
num_metal = 0
coor_mol = np.zeros((num_mol,2))
coor_metal = np.zeros((num_metal,2))
################### Distribute the molecules ####################
for i in range(0,num_mol):
    state = True
    while state == True:
        ind_x = rd.randint(0, latt_len-1)
        ind_y = rd.randint(0, latt_len-1)
        pos_current = get_coor_mol([ind_x,ind_y])
        if is_occupied(pos_current, lattice) == False:
            set_element([ind_x,ind_y],1,i,coor_mol,lattice)
            set_element([ind_x,ind_y],i,i,coor_mol,lattice_num)
            state = False
#print lattice_num
print "molecules are distributed, %s" % time.clock()

for i in range(0,num_metal):
    state = True
    while state == True:
        ind_x = rd.randint(0, latt_len-1)
        ind_y = rd.randint(0, latt_len-1)
        if is_occupied([ind_x,ind_y], lattice) == False:
            set_element([ind_x,ind_y],1,i,coor_metal,lattice)
            set_element([ind_x,ind_y],latt_len,i,coor_metal,lattice_num)
            state = False

            
################### DO THE MONTE CARLO SIMULATION ############################
Sim_Enabled =  False
print "Simulation begins..."
if Sim_Enabled == True:
    count = 0
    while count < TOTAL_RUN:
        ind_element = rd.randint(0, num_mol+num_metal-1)
        if ind_element < num_mol:
            energy_current = cal_energy_mol(coor_mol[ind_element,:],coor_mol,lattice_num)
            #### It is optional whether to remove the selected molecules before go on ######
            state = True
            while state == True:
                new_mol_pos = [rd.randint(0,latt_len-1), rd.randint(0,latt_len-1)]
                if is_occupied(new_mol_pos, lattice) == False: #and is_forbidden(new_mol_pos,lattice) == False:
                    energy_new = cal_energy_mol(new_mol_pos,coor_mol,lattice_num)
                    p = min(math.exp(-(energy_new - energy_current)),1)
                    if p > rd.random():
                        set_mol(coor_mol[ind_element,:],0,ind_element,coor_mol,lattice)
                        set_mol(coor_mol[ind_element,:],0,ind_element,coor_mol,lattice_num)
                        set_mol(new_mol_pos,1,ind_element,coor_mol,lattice)
                        set_mol(new_mol_pos,ind_element,ind_element,coor_mol,lattice_num)
                state = False
        else:
            ind_element = ind_element - num_metal
            energy_current = cal_energy_metal(coor_metal[ind_element,:],coor_metal,lattice)
            state = True
            while state == True:
                new_metal_pos = [rd.randint(0,latt_len-1), rd.randint(0,latt_len-1)]
                if is_occupied(new_mol_pos, lattice) == False:
                    energy_new = cal_energy_metal(new_metal_pos,coor_metal,lattice)
                    p = min(math.exp(-(energy_new - energy_current)),1)
                    if p > rd.random():
                        pass



        count = count + 1
        if count%(TOTAL_RUN/10) == 0:
            print "number of run: %d / 10, costed time: %f" % (count/(TOTAL_RUN/10), time.clock())





plt.figure(1)
plt.imshow(lattice)
plt.figure(2)
plt.imshow(lattice_num)
plt.show()





