# Monte Carlo

import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import copy
import time
import Image

# global lattice_size
LATTICE_CONSTANT = 1
TOTAL_RUN = 100000000







def sir(x,range):
	#global latt_len
	if x > range-1:
		x = x - range
	elif x < 0:
		x = x + range
	return x


def get_coor_mol(input_coor):
	#if input_coor.size == 2:
		#x = input_coor[0,0]
		#y = input_coor[0,1]
	#else:
	x = input_coor[0]
	y = input_coor[1]
	coor = np.zeros((5,2))
	coor[0,:] = [x, y]
	coor[1,:] = [sir(x-1,latt_len), y]
	coor[2,:] = [x, sir(y+1,latt_len)]
	coor[3,:] = [sir(x+1,latt_len), y]
	coor[4,:] = [x, sir(y-1,latt_len)]
	return coor

def set_element(input_coor,op,id_ele,coor,latt):
    ele_length = input_coor.size
    if ele_length == 0:
        print "Input is invalid!"
        return
    elif ele_length == 2:
        latt[input_coor[0],input_coor[1]] = op
        coor[id_ele,:] = input_coor
        #print input_coor
    elif ele_length > 2:
        for i in range(0,ele_length/2):
            latt[input_coor[i,0],input_coor[i,1]] = op
        coor[id_ele,:] = input_coor[0,:]

def is_occupied(input_coor,latt):
	#print input_coor
	ele_length = input_coor.size
	if ele_length == 2:
		if latt[input_coor[0],input_coor[1]] != 0:
			return True
	elif input_coor.shape[0] > 2:
		for i in range(0,ele_length/2):
			if latt[input_coor[i,0],input_coor[i,1]] != 0:
				return True
	return False

def is_forbidden(input_coor,latt):
	#return False
	direct = np.zeros((4,2))
	direct[0,:] = [-1,0]
	direct[1,:] = [0,+1]
	direct[2,:] = [+1,0]
	direct[3,:] = [0,-1]
	ele_length = input_coor.size
	if ele_length == 2:
		count = []
		#print input_coor
		pos_around = get_coor_mol(input_coor)[1:,:]
		pos_around2 = pos_around + direct
		for i in range(0,4):
			if latt[pos_around[i,0],pos_around[i,1]] == latt_len:
				return True
			elif latt[pos_around[i,0],pos_around[i,1]] != 0 and latt[pos_around[i,0],pos_around[i,1]] != latt_len: 
				if latt[pos_around[i,0],pos_around[i,1]] == \
						latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)]:
							count.append(i)
		#print count
		if len(count) > 2:
			#print count
			return True
		elif len(count) == 2:
			if count[0] == 0 and count[1] == 2:
				return False
			elif count[0] == 1 and count[1] == 3:
				return False
			else:
				#print count
				return True
	elif ele_length > 2:
		#print input_coor
		pos_around = input_coor[1:,:]
		pos_around2 = pos_around + direct
		for i in range(0,4):
			if latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)] == latt_len:
				plus1 = pos_around2[i,:] + direct[sir(i+1,4),:]
				plus2 = plus1 + direct[sir(i+1,4),:]
				#print "test is forbiden"
				#print pos_around
				#print pos_around2
				#print plus1
				#print plus2
				minus1 = pos_around2[i,:] + direct[sir(i-1,4),:]
				minus2 = minus1 + direct[sir(i-1,4),:]
				if latt[sir(plus1[0],latt_len),sir(plus1[1],latt_len)] == latt[sir(plus2[0],latt_len),sir(plus2[1],latt_len)] \
						or latt[sir(minus1[0],latt_len),sir(minus1[1],latt_len)] == latt[sir(minus2[0],latt_len),sir(minus2[1],latt_len)]:
							return True
	return False

def cal_energy_mol(coor,coor_mol,cenergy,venergy,latt):
    # First detect if there are molecules around
    energy = 0
    temp = get_coor_mol(coor)
    direct = np.zeros((4,2))
    direct[0,:] = [-1,0]
    direct[1,:] = [0,+1]
    direct[2,:] = [+1,0]
    direct[3,:] = [0,-1]
    # Get the points around this molecule
    #pos_around = np.zeros((4,2))
    pos_around = temp[1:5,:] + direct
    pos_around2 = pos_around + direct
    for i in range(0,4):
        if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] == latt_len:
            energy = energy - cenergy
        elif latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != 0:
            if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != \
                    latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)]: 
			energy = energy - venergy
    return energy

def cal_energy_metal(coor,coor_metal,mc_energy,latt):
    energy = 0
    direct = np.zeros((4,2))
    direct[0,:] = [-1,0]
    direct[1,:] = [0,+1]
    direct[2,:] = [+1,0]
    direct[3,:] = [0,-1]
    pos_around = get_coor_mol(coor)[1:,:]
    pos_around2 = pos_around + direct
    #print pos_around
    for i in range(0,4):
        if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != 0 and latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != latt_len: 
            if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] == \
                    latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)]:
                        energy = energy - mcenergy 
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

def experiment(total_run, latt_len, num_mol, num_metal, cenergy, venergy, mcenergy):
	print "Begin, %s" % time.clock()
	#print "Total run: %d" % total_run
	lattice = np.zeros((latt_len,latt_len))
	lattice_num = np.zeros((latt_len,latt_len))
	
	coor_mol = np.zeros((num_mol,2))
	coor_metal = np.zeros((num_metal,2))
	################### Distribute the molecules ####################
	#print "Distributing molecules..."
	for i in range(0,num_mol):
	    state = True
	    while state == True:
	        ind_x = rd.randint(0, latt_len-1)
	        ind_y = rd.randint(0, latt_len-1)
	        pos_current = get_coor_mol(np.array([ind_x,ind_y]))
	        if is_occupied(pos_current, lattice) == False and is_forbidden(pos_current,lattice_num) == False:
	            set_element(pos_current,2,i,coor_mol,lattice)
	            set_element(pos_current,i,i,coor_mol,lattice_num)
	            state = False
	#print "Molecules are distributed, %s" % time.clock()
	
	#print "Distributing metals..."
	for i in range(0,num_metal):
	    state = True
	    while state == True:
	        ind_x = rd.randint(0, latt_len-1)
	        ind_y = rd.randint(0, latt_len-1)
	        pos_current = np.array([ind_x,ind_y])
	        if is_occupied(pos_current, lattice) == False and is_forbidden(pos_current, lattice_num) == False:
	            set_element(pos_current,1,i,coor_metal,lattice)
	            set_element(pos_current,latt_len,i,coor_metal,lattice_num)
	            state = False
	#print "Metals are distributed..."
	            
	################### DO THE MONTE CARLO SIMULATION ############################
	#print "Simulation begins..."
	count = 0
	while count < total_run:
		ind_element = rd.randint(0, num_mol+num_metal-1)
		#ind_element = rd.randint(0, num_mol-1)
		if ind_element < num_mol:
			energy_current = cal_energy_mol(coor_mol[ind_element,:],coor_mol,cenergy,venergy,lattice_num)
			#### It is optional whether to remove the selected molecules before go on ######
			state = True
			while state == True:
				new_mol_pos = np.array([rd.randint(0,latt_len-1), rd.randint(0,latt_len-1)])
				pos_new = get_coor_mol(new_mol_pos)
				if is_occupied(pos_new, lattice) == False and is_forbidden(pos_new,lattice_num) == False:
					pos_old = get_coor_mol(coor_mol[ind_element,:])
					energy_new = cal_energy_mol(new_mol_pos,coor_mol,cenergy,venergy,lattice_num)
					p = min(math.exp(-(energy_new - energy_current)),1)
					if p > rd.random():
						set_element(pos_old,0,ind_element,coor_mol,lattice)
						set_element(pos_old,0,ind_element,coor_mol,lattice_num)
						set_element(pos_new,2,ind_element,coor_mol,lattice)
						set_element(pos_new,ind_element,ind_element,coor_mol,lattice_num)
				state = False
		else:
			ind_element = ind_element - num_mol
			#print ind_element
			energy_current = cal_energy_metal(coor_metal[ind_element,:],coor_metal,mcenergy,lattice_num)
			state = True
			while state == True:
				new_metal_pos = np.array([rd.randint(0,latt_len-1), rd.randint(0,latt_len-1)])
				if is_occupied(new_metal_pos, lattice) == False and is_forbidden(new_metal_pos,lattice_num) == False:
					old_metal_pos = coor_metal[ind_element,:]
					energy_new = cal_energy_metal(new_metal_pos,coor_metal,mcenergy,lattice_num)
					p = min(math.exp(-(energy_new - energy_current)),1)
					if p > rd.random():
						set_element(old_metal_pos,0,ind_element,coor_metal,lattice)
						set_element(old_metal_pos,0,ind_element,coor_metal,lattice_num)
						set_element(new_metal_pos,1,ind_element,coor_metal,lattice)
						set_element(new_metal_pos,latt_len,ind_element,coor_metal,lattice_num)
					state = False
		count = count + 1
		#if count%(total_run/10) == 0:
		#	print "number of run: %d / 10, costed time: %f" % (count/(total_run/10), time.clock())
	#lattice_num = lattice_num.astype(np.uint8)
	#im = Image.fromarray(lattice_num)
	#im.save("results/%d-%d-%d-%d-%d.jpeg" % (num_mol,num_metal,cenergy,venergy,mcenergy))
	plt.imsave("results1\%d-%d-%d-%d-%d-%d-%.1e.png" % (num_mol,num_metal,cenergy,venergy,mcenergy,latt_len,total_run),lattice,[0,2])
	np.savetxt('results1\%d-%d-%d-%d-%d-%d-%.1e.txt' % (num_mol,num_metal,cenergy,venergy,mcenergy,latt_len,total_run), lattice_num, fmt='%i', delimiter=',', comments = '(%.1e-%d-%d-%d-%d-%d' % (total_run,num_mol,num_metal \
			,cenergy, venergy, mcenergy))
	print "End, %s" % time.clock()

total_run = 10000

latt_len = 50
num_mol = 50
num_metal = 10

cenergy = 25
venergy = 1
mcenergy = 25
while num_metal < 150:
	while venergy < 10:
		 experiment(total_run,latt_len,num_mol,num_metal,cenergy,venergy,mcenergy)
		 venergy = venergy + 1
	num_metal = num_metal + 10

#print "Simulation is done! costed time: %f" % (time.clock())

#plt.figure(1)
#plt.imshow(lattice)

#plt.figure(2)
#plt.imshow(lattice_num)
#plt.show()

#np.savetxt('\results\%d-%d.txt' % , lattice, fmt='%i', delimiter=',')






