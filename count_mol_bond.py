# Python file to calculate the bonding of the system
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy



# load the data
dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results14"
os.chdir(dname)
latt_len = 80
files = os.listdir(dname)
line = files[47]
print line
lattice = np.loadtxt(line, delimiter=',',skiprows=1)
##################################################
def prepro(lattice):
	latt_len = lattice.shape[0]
	mol = np.transpose(np.array(np.where(lattice == 3)))
	output = np.zeros((latt_len,latt_len))
	#print mol.shape[0]
	for i in range(mol.shape[0]):
		output[mol[i][0],mol[i][1]] = i + 1
	metal = np.where(lattice == 1)
	output[metal] = 1000
	return output

def sir(x,range):
	#global latt_len
	if x > range-1:
		x = x - range
	elif x < 0:
		x = x + range
	return x

def find_neighbour(coor, direct, lattice):
	
	#print "input coordinate for find_neighbour: " + str(coor)
	neighbours = []
	latt_len = lattice.shape[0]
	x = coor[0]
	y = coor[1]
	for i in range(0,direct.shape[0]):
		temp_coor1 = [sir(x+direct[i,0],latt_len),sir(y+direct[i,1],latt_len)]
		if lattice[temp_coor1[0]][temp_coor1[1]] == 1000:
			temp_coor2 = [sir(temp_coor1[0]+direct[i,0],latt_len),sir(temp_coor1[1]+direct[i,1],latt_len)]
			if lattice[temp_coor2[0]][temp_coor2[1]] != 0 and lattice[temp_coor2[0]][temp_coor2[1]] != 1000:
				neighbours.append(lattice[temp_coor2[0],temp_coor2[1]]-1)
	return neighbours
				
def cluster(lattice):
	direct = np.zeros((4,2))
	direct[0,:] = [-2,0]
	direct[1,:] = [0,+2]
	direct[2,:] = [+2,0]
	direct[3,:] = [0,-2]
	latt = deepcopy(lattice)
	latt = prepro(latt)
	mol = np.transpose(np.array(np.where(lattice == 3)))
	num_mol = mol.shape[0]
	output = np.array(range(1,num_mol+1))
	output_cluster = []
	if output.shape[0] != num_mol:
		print "WARNING!!! WRONG INITIALIZATION OF ELEMENT!!!"
	for i in range(0,num_mol):
		if output[i] == i+1:
			temp = [i] # list that store the connected elements
			templist = [i] # Queue, first in first out
			while len(templist) != 0:
				# pop out the first element
				current_index = templist.pop(0)
				# find out the neighbours of this element
				connected_neighbours = find_neighbour(mol[current_index,:],direct,latt)
				#print "connected_neighbours " + str(connected_neighbours)
				# loop through all the neighbours, if it is not the first element and if it has been already asigned
				for j in range(len(connected_neighbours)):
					if connected_neighbours[j] != i and output[connected_neighbours[j]] == connected_neighbours[j]+1:
						# change the label of these elements to the current_index
						output[connected_neighbours[j]] = i+1
						# add it to the check list
						templist.append(connected_neighbours[j])
						# add it to the temp list contains all the elements of this cluster
						temp.append(int(connected_neighbours[j]))
						# change the label of these elements to the current_index
						latt[mol[connected_neighbours[j]][0]][mol[connected_neighbours[j]][1]] = i+1
						
			output_cluster.append(temp)
	return (output_cluster,latt)

def corr_num(coor, direction, lattice):
	x = coor[0]
	y = coor[1]
	latt_len = lattice.shape[0]
	cal_num = direction.shape[0]
	#print "cal_num " + str(cal_num)
	count = 0
	for i in range(cal_num):
		temp = lattice[sir(x+direction[i][0],latt_len)][sir(y+direction[i][1],latt_len)]
		#print "temp " + str(temp)
		#if temp != 0 and temp != 1000:
		if temp == 3:
			count = count + 1
	return count
			

# auto-correlation for different mode: 0 for dense-packed, 1 for 1D, 2 for 2D		
def auto_correlate(mode, threshold, element, index, lattice):
	latt = deepcopy(lattice)
	latt = prepro(latt)
	# 1. cluster
	#output,latt = cluster(lattice)
	# 2. auto-correlation
	th = threshold # elements with auto-correlation larger than threshold are chosen
	mol_count = 0
	# store the picked elements for debugging
	mol_list = []
	# get the positions of elements
	#mol = np.transpose(np.array(np.where(lattice == 3)))
	mol = element
	num_mol = mol.shape[0]
	#print num_mol
	if mode == 2:
		direction = np.zeros((8,2))
		direction[0,:] = [-4,0]
		direction[1,:] = [-4,+4]
		direction[2,:] = [0,+4]
		direction[3,:] = [+4,+4]
		direction[4,:] = [+4,0]
		direction[5,:] = [+4,-4]
		direction[6,:] = [0,-4]
		direction[7,:] = [-4,-4]
		# go through each cluster
		for i in range(num_mol):
			temp_count = corr_num(mol[i,:],direction,lattice)
			#print "entry " + str(i) + "  temp_count " + str(temp_count)
			if temp_count > th:
				mol_count = mol_count + 1
				mol_list.append(index[i])
				latt[mol[i][0]][mol[i][1]] = 1001
	return (mol_count,mol_list)
			
plt.figure()	
plt.imshow(lattice)
plt.show()
output,latt = cluster(lattice)
for i in range(len(output)):
	print "new: " + str(len(output[i])) + "  " + str(output[i])

mol = np.transpose(np.array(np.where(lattice == 3)))
num_mol = mol.shape[0]
new_mol = []
for i in range(len(output)):
	if len(output[i]) > 3:
		ind = np.transpose(np.array(output[i]))
		#print ind
		
		ind_num = ind.shape[0]
		#print ind_num
		ele = mol[ind,:]
		count, mols = auto_correlate(2, 2, ele, ind, lattice)
	new_mol.append(mols)
print "###############################"
for i in range(len(new_mol)):
	print "new: " + str(len(new_mol[i])) + "  " + str(new_mol[i])


new_latt = np.zeros((lattice.shape[0],lattice.shape[0]))
for i in range(len(new_mol)):
	for j in range(len(new_mol[i])):
		#print new_mol[i][j]
		new_latt[mol[new_mol[i][j]][0]][mol[new_mol[i][j]][1]] = (i+1)*3
	
legs = np.where(lattice == 2)
lattice[legs] = 0

metal = np.where(lattice == 1)
lattice[metal] = 0


plt.imshow(new_latt)
plt.show()








#for line in files:
	#if (line[-4:] == '.txt'):
		#namedata = line[0:-4].strip().split('-')
		#num_metal = int(namedata[3])
		#num_metal_ind = num_metal/50 - 1
		#cenergy = int(namedata[4])
		#cenergy_ind = (cenergy -1)/2
		#f = open(line, 'r')
		#headdata = f.readline().strip().split(',')
		#f.close()
		#cbondvalue[cenergy_ind][num_metal_ind] = float(headdata[0])#/float(num_metal)
		#vbondvalue[cenergy_ind][num_metal_ind] = int(headdata[1])
		#energyvalue[cenergy_ind][num_metal_ind] = int(headdata[2])

