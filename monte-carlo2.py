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
lattice_size = 40
LATTICE_CONSTANT = 1
TOTAL_RUN = 600000000
print "Total run: %d" % TOTAL_RUN 
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
def det_neighbour(current_x, current_y, angle):
	global lattice_size
	# angle: 0,1,2,3,4,5
	# angle = 0 means 0 degree, as angle increases, count it anti-clockwise
	i = current_x
	j = current_y
	if angle == 0:
		temp = [i,j+1]
	elif angle == 1:
		temp = [i-1,j+1]
	elif angle == 2:
		temp = [i-1,j]
	elif angle == 3:
		temp = [i,j-1]
	elif angle == 4:
		temp = [i+1,j-1]
	else:
		temp = [i+1,j]
	if temp[0] < 0:
		temp[0] = temp[0] + lattice_size
	elif temp[0] > lattice_size-1:
		temp[0] = temp[0] - lattice_size
	if temp[1] < 0:
		temp[1] = temp[1] + lattice_size
	elif temp[1] > lattice_size-1:
		temp[1] = temp[1] - lattice_size
	return temp

def keep_within_5(var):
	#return_value
	if var > 5:
		return_value = var - 6
	elif var < 0:
		return_value = var + 6
	else:
		return_value = var
	return return_value

def det_anti(angle):
	if angle > 2:
		angle_anti = angle - 3
	else:
		angle_anti = angle + 3
	return angle_anti
		
def is_occupied(lat, co):

	ele_length = len(co)
	for i in range(0,ele_length):
		if lat[co[i][0],co[i][1]] == 0:
			continue
		else:
			return True
	return False

def is_forbidden(lat, element):
	ele_length = len(element)
	if ele_length == 2:
		# it is a metal atom
		# check if there is already a molecule around
		ind_list =[]
		for i in range(0,5):
			temp1 = det_neighbour(element[0],element[1],i)
			temp2 = det_neighbour(temp1[0],temp1[1],i)
			if lattice[temp1[0],temp1[1]] == 2 and lattice[temp2[0],temp2[1]] == 3:
				ind_list.append(i)
		int_count = len(ind_list)
		if int_count <= 1:
			return False
		elif int_count > 3:
			return True
		elif int_count == 2:
			if abs(ind_list[0]-ind_list[1]) == 1 or abs(ind_list[0]-ind_list[1]) == 5:
				return True
			else:
				return False
		else:
			if (ind_list[2]-ind_list[1]) == 2 and (ind_list[1]-ind_list[0]) == 2:
				return False
			else:
				return True
	else:
		# it is molecule
		ind_mo1 = det_neighbour(element[0],element[1],element[2])
		# first possible metal position
		temp1 = det_neighbour(ind_mo1[0],ind_mo1[1],element[2])
		# if there is a metal atom
		if lat[temp1[0],temp1[1]] == 1:
			anti_angle = det_anti(element[2])
			anti_angle_nei1 = keep_within_5(anti_angle - 1)
			ind_mo1_nei1 = det_neighbour(temp1[0],temp1[1],anti_angle_nei1)
			ind_mo1_nei11 = det_neighbour(ind_mo1_nei1[0],ind_mo1_nei1[1],anti_angle_nei1)
			if lat[ind_mo1_nei1[0],ind_mo1_nei1[1]] == 2 and lat[ind_mo1_nei11[0],ind_mo1_nei11[1]] == 3:
				return True
			
			anti_angle_nei2 = keep_within_5(anti_angle + 1)
			ind_mo1_nei2 = det_neighbour(temp1[0],temp1[1],anti_angle_nei2)
			ind_mo1_nei22 = det_neighbour(ind_mo1_nei2[0],ind_mo1_nei2[1],anti_angle_nei2)
			if lat[ind_mo1_nei2[0],ind_mo1_nei2[1]] == 2 and lat[ind_mo1_nei22[0],ind_mo1_nei22[1]] == 3:
				return True
				
		ind_mo2 = det_neighbour(element[0],element[1],element[3])
		# second possible metal position
		temp2 = det_neighbour(ind_mo2[0],ind_mo2[1],element[3])
		if lat[temp2[0],temp2[1]] == 1:
			anti_angle = det_anti(element[3])
			anti_angle_nei1 = keep_within_5(anti_angle - 1)
			ind_mo2_nei1 = det_neighbour(temp2[0],temp2[1],anti_angle_nei1)
			ind_mo2_nei11 = det_neighbour(ind_mo2_nei1[0],ind_mo2_nei1[1],anti_angle_nei1)
			if lat[ind_mo2_nei1[0],ind_mo2_nei1[1]] == 2 and lat[ind_mo2_nei11[0],ind_mo2_nei11[1]] == 3:
				return True
			
			anti_angle_nei2 = keep_within_5(anti_angle + 1)
			ind_mo2_nei2 = det_neighbour(temp2[0],temp2[1],anti_angle_nei2)
			ind_mo2_nei22 = det_neighbour(ind_mo2_nei2[0],ind_mo2_nei2[1],anti_angle_nei2)
			if lat[ind_mo2_nei2[0],ind_mo2_nei2[1]] == 2 and lat[ind_mo2_nei22[0],ind_mo2_nei22[1]] == 3:
				return True
		return False
			
			
def cal_energy_mol(ind_mo, lattice):
	energy = 0
	ind_mo1 = det_neighbour(ind_mo[0],ind_mo[1],ind_mo[2])
	temp1 = det_neighbour(ind_mo1[0],ind_mo1[1],ind_mo[2])
	if lattice[temp1[0],temp1[1]] == 1:
		energy = energy - 10
	ind_mo2 = det_neighbour(ind_mo[0],ind_mo[1],ind_mo[3])
	temp2 = det_neighbour(ind_mo2[0],ind_mo2[1],ind_mo[3])
	if lattice[temp2[0],temp2[1]] == 1:
		energy = energy - 10
	return energy

def cal_energy_metal(ind_me, lattice):
	energy = 0
	for i in range(0,5):
		temp1 = det_neighbour(ind_me[0],ind_me[1],i)
		temp2 = det_neighbour(temp1[0],temp1[1],i)
		if lattice[temp1[0],temp1[1]] == 2 and lattice[temp2[0],temp2[1]] == 3:
			energy = energy - 10
	return energy


print "clock1:%s" % time.clock() 
################################### define the lattice
lattice = np.zeros((lattice_size, lattice_size))
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


number_molecule = 72 # represent by 2,3,2
number_metal = 48 # represent by 1

angle_group = np.array([])

############################## put all the molecules on the lattice
molecules = np.zeros((number_molecule,4))
for i in range(0,number_molecule):
	state = True
	while state == True:
		ind_x = rd.randint(0, lattice_size-1)
		ind_y = rd.randint(0, lattice_size-1)
		angle = rd.randint(0,5)
		angle_anti = det_anti(angle)
		points_t1 = det_neighbour(ind_x,ind_y,angle)
		points_t2 = det_neighbour(ind_x,ind_y,angle_anti)
		points_temp = [[ind_x,ind_y],points_t1,points_t2]
		#print points_temp
		if is_occupied(lattice, points_temp) == False:
			lattice[points_temp[0][0],points_temp[0][1]] = 3
			lattice[points_temp[1][0],points_temp[1][1]] = 2
			lattice[points_temp[2][0],points_temp[2][1]] = 2
			state = False
			#molecules[i,:] = [ind_x,ind_y,points_t1[0],points_t1[1],points_t2[0],points_t2[1]]
			molecules[i,:] = [ind_x, ind_y, angle, angle_anti]
	
################################ put all the metals on the lattice
metals = np.zeros((number_metal,2))
for i in range(0,number_metal):
	state = True
	while state == True:
		ind_x = rd.randint(0, lattice_size-1)
		ind_y = rd.randint(0, lattice_size-1)
		temp_points = [[ind_x,ind_y]]
		if is_occupied(lattice, temp_points) == False and is_forbidden(lattice,[ind_x,ind_y]) == False:
			lattice[ind_x,ind_y] = 1
			metals[i,:] = [ind_x, ind_y]
			state = False
			

print "clock2:%s" % time.clock()
temp1, temp2 = np.where(lattice == 3)
print len(temp1)
lattice2 = copy.deepcopy(lattice)
np.savetxt("text.txt", lattice, fmt="%d", delimiter=",")
############################### begin the main loop
i = 0
while i < TOTAL_RUN:
	# randomly pick a molecule and a metal
	ind_element = rd.randint(0, number_molecule + number_metal-1)
	# if choose a molecule
	if ind_element < number_molecule:
		#lattice_temp = copy.deepcopy(lattice)
		# calculate the energy of current configuration
		energy_current = cal_energy_mol(molecules[ind_element,:],lattice)
		# calculate the current molecule positions
		points_current = [molecules[ind_element,0:2].tolist(),\
			det_neighbour(molecules[ind_element,0],molecules[ind_element,1],molecules[ind_element,2]),\
			det_neighbour(molecules[ind_element,0],molecules[ind_element,1],molecules[ind_element,3])]
		#print points_current
		# assign the updated value to the temp lattice
		lattice[points_current[0][0],points_current[0][1]] = 0
		lattice[points_current[1][0],points_current[1][1]] = 0
		lattice[points_current[2][0],points_current[2][1]] = 0
		state = True
		while state == True:
			new_mol_pos = [rd.randint(0, lattice_size-1),rd.randint(0, lattice_size-1)]
			new_mol_angle = rd.randint(0,5)
			new_mol_anti_angle = det_anti(new_mol_angle)
			new_mol_points = [new_mol_pos[0], new_mol_pos[1], new_mol_angle, new_mol_anti_angle]
			points_temp = [new_mol_pos, det_neighbour(new_mol_pos[0],new_mol_pos[1],new_mol_angle),\
				det_neighbour(new_mol_pos[0],new_mol_pos[1],new_mol_anti_angle)]
			#print points_temp
			if is_occupied(lattice, points_temp) == False and is_forbidden(lattice, new_mol_points) == False:
				# further updated the temp lattice
				energy_new = cal_energy_mol(new_mol_points,lattice)
				p = min(math.exp(-(energy_new - energy_current)),1)
				if p > rd.random():
					lattice[points_temp[0][0],points_temp[0][1]] = 3
					lattice[points_temp[1][0],points_temp[1][1]] = 2
					lattice[points_temp[2][0],points_temp[2][1]] = 2
					#lattice = copy.deepcopy(lattice_temp)
					#lattice = lattice_temp
					molecules[ind_element,:] = [new_mol_pos[0],new_mol_pos[1], new_mol_angle, new_mol_anti_angle]
				else:
					lattice[points_current[0][0],points_current[0][1]] = 3
					lattice[points_current[1][0],points_current[1][1]] = 2
					lattice[points_current[2][0],points_current[2][1]] = 2
				state = False
				
	# if choose a metal
	#if ind_element >= number_molecule:
	else:	
		ind_element = ind_element - number_molecule
		#print ind_element
		#lattice_temp = copy.deepcopy(lattice)
		# calculate the energy of current configuration
		energy_current = cal_energy_metal(metals[ind_element,:],lattice)
		# calculate the current molecule positions
		points_current = metals[ind_element,:].tolist()
		#print points_current
		# assign the updated value to the temp lattice
		lattice[points_current[0],points_current[1]] = 0
		state = True
		while state == True:
			new_metal_pos = [rd.randint(0, lattice_size-1),rd.randint(0, lattice_size-1)]
			points_temp = [new_metal_pos]
			#print points_temp
			if is_occupied(lattice, points_temp) == False and is_forbidden(lattice, new_metal_pos) == False:
				# further updated the temp lattice
				energy_new = cal_energy_metal(new_metal_pos,lattice)
				p = min(math.exp(-(energy_new - energy_current)),1)
				if p > rd.random():
					lattice[points_temp[0][0],points_temp[0][1]] = 1
					#lattice = copy.deepcopy(lattice_temp)
					#lattice = lattice_temp
					metals[ind_element,:] = [new_metal_pos[0],new_metal_pos[1]]
				else:
					lattice[points_current[0],points_current[1]] = 1
				state = False
	i= i + 1
	if i%(TOTAL_RUN/10) == 0:
		print "number of run: %d / 10, costed time : %f" % (i/(TOTAL_RUN/10), time.clock())
print "clock3:%s" % time.clock()
######################### display
#plt.figure()
##plt.subplot(1, 2, 1)
#for i in range(0, lattice_size):
	#plt.scatter(lattx[i,:], latty[i,:])
#temp1, temp2 = np.where(lattice2 == 3)
#temp3, temp4 = np.where(lattice2 == 2)
#temp5, temp6 = np.where(lattice2 == 1)
#print len(temp1)
#plt.scatter(lattx[temp1,temp2],latty[temp1,temp2],s = 70,c = "r",marker = "s")
#plt.scatter(lattx[temp3,temp4],latty[temp3,temp4],s = 60,c = "y",marker = "s")
#plt.scatter(lattx[temp5,temp6],latty[temp5,temp6],s = 60,c = "k",marker = "s")
##im = plt.imshow(lattice)
#plt.show()
#plt.subplot(1, 2, 2)
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
