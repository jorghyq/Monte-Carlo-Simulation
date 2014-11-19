# Python file to calculate the bonding of the system
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

#direct0 = np.zeros((8,2))
#direct0[0,:] = [-2,+1]
#direct0[1,:] = [-1,+2]
#direct0[2,:] = [+1,+2]
#direct0[3,:] = [+2,+1]
#direct0[4,:] = [+2,-1]
#direct0[5,:] = [+1,-2]
#direct0[6,:] = [-1,-2]
#direct0[7,:] = [-2,-1]
#direct1 = np.zeros((12,2))
#direct1[0,:] = [-2,+1]
#direct1[1,:] = [-1,+2]
#direct1[2,:] = [+1,+2]
#direct1[3,:] = [+2,+1]
#direct1[4,:] = [+2,-1]
#direct1[5,:] = [+1,-2]
#direct1[6,:] = [-1,-2]
#direct1[7,:] = [-2,-1]
#direct1[8,:] = [-2,0]
#direct1[9,:] = [0,+2]
#direct1[10,:] = [+2,0]
#direct1[11,:] = [0,-2]
#direct2 = np.zeros((4,2))
#direct2[0,:] = [-2,0]
#direct2[1,:] = [0,+2]
#direct2[2,:] = [+2,0]
#direct2[3,:] = [0,-2]

# load the data

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

def find_neighbour(mode, coor, direct, lattice):
	neighbours = []
	latt_len = lattice.shape[0]
	x = coor[0]
	y = coor[1]
	if mode == 2:
		
		for i in range(0,direct.shape[0]):
			temp_coor1 = [sir(x+direct[i,0],latt_len),sir(y+direct[i,1],latt_len)]
			if lattice[temp_coor1[0]][temp_coor1[1]] == 1000:
				temp_coor2 = [sir(temp_coor1[0]+direct[i,0],latt_len),sir(temp_coor1[1]+direct[i,1],latt_len)]
				if lattice[temp_coor2[0]][temp_coor2[1]] != 0 and lattice[temp_coor2[0]][temp_coor2[1]] != 1000:
					neighbours.append(lattice[temp_coor2[0],temp_coor2[1]]-1)
	elif mode == 0:
		for i in range(0,direct.shape[0]):
			temp_coor1 = [sir(x+direct[i,0],latt_len),sir(y+direct[i,1],latt_len)]
			if lattice[temp_coor1[0]][temp_coor1[1]] != 0 and lattice[temp_coor1[0]][temp_coor1[1]] != 1000:
				neighbours.append(lattice[temp_coor1[0],temp_coor1[1]]-1)
	elif mode == 1:
		
		direct_dense = direct[0:8,:]
		direct_2d = direct[8:12,:]
		for i in range(0,direct_dense.shape[0]):
			temp_coor1 = [sir(x+direct_dense[i,0],latt_len),sir(y+direct_dense[i,1],latt_len)]
			if lattice[temp_coor1[0]][temp_coor1[1]] != 0 and lattice[temp_coor1[0]][temp_coor1[1]] != 1000:
				if len(find_neighbour(0,temp_coor1,direct_dense,lattice)) < 3:
					neighbours.append(lattice[temp_coor1[0],temp_coor1[1]]-1)
		for i in range(0,direct_2d.shape[0]):
			temp_coor1 = [sir(x+direct_2d[i,0],latt_len),sir(y+direct_2d[i,1],latt_len)]
			if lattice[temp_coor1[0]][temp_coor1[1]] == 1000:
				temp_coor2 = [sir(temp_coor1[0]+direct_2d[i,0],latt_len),sir(temp_coor1[1]+direct_2d[i,1],latt_len)]
				if lattice[temp_coor2[0]][temp_coor2[1]] != 0 and lattice[temp_coor2[0]][temp_coor2[1]] != 1000:
					if len(find_neighbour(2,temp_coor2,direct_2d,lattice)) < 3:
						neighbours.append(lattice[temp_coor2[0],temp_coor2[1]]-1)					
	return neighbours

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
	# 1. cluster
	#output,latt = cluster(lattice)
	# 2. auto-correlation
	latt_len = lattice.shape[0]
	latt = np.zeros((latt_len,latt_len))
	for i in range(element.shape[0]):
		latt[tuple(element[i,:])] = 3
	#metal = np.where(lattice == 1)
	#latt[metal] = 1
	th = threshold # elements with auto-correlation larger than threshold are chosen
	mol_count = 0
	# store the picked elements for debugging
	mol_list = []
	# get the positions of elements
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
			temp_count = corr_num(mol[i,:],direction,latt)

			#print "entry " + str(i) + "  temp_count " + str(temp_count)
			if temp_count > th:
				mol_count = mol_count + 1
				mol_list.append(index[i])
				#latt[mol[i][0]][mol[i][1]] = 1001
	elif mode == 0:
		direct0 = np.zeros((16,2))
		direct0[0,:] = [-2,+1]
		direct0[1,:] = [-1,+2]
		direct0[2,:] = [+1,+2]
		direct0[3,:] = [+2,+1]
		direct0[4,:] = [+2,-1]
		direct0[5,:] = [+1,-2]
		direct0[6,:] = [-1,-2]
		direct0[7,:] = [-2,-1]
		# consider the neighbours in the diagonal
		direct0[8,:] = [+3,-1]
		direct0[9,:] = [+3,+1]
		direct0[10,:] = [-3,+1]
		direct0[11,:] = [-3,-1]
		direct0[12,:] = [+1,-3]
		direct0[13,:] = [+1,+3]
		direct0[14,:] = [-1,+3]
		direct0[15,:] = [-1,-3]
		for i in range(num_mol):
			temp_count = corr_num(mol[i,:],direct0,latt)
			#if index[i] == 182:
			#	print "temp_count 108 : " + str(temp_count)
			#print "entry " + str(i) + "  temp_count " + str(temp_count)
			if temp_count > th:
				mol_count = mol_count + 1
				mol_list.append(index[i])
	elif mode == 1:
		direct1 = np.zeros((28,2))
		direct1[0,:] = [-2,+1]
		direct1[1,:] = [-1,+2]
		direct1[2,:] = [+1,+2]
		direct1[3,:] = [+2,+1]
		direct1[4,:] = [+2,-1]
		direct1[5,:] = [+1,-2]
		direct1[6,:] = [-1,-2]
		direct1[7,:] = [-2,-1]
		# consider the neighbours in the diagonal
		direct1[8,:] = [-4,0]
		direct1[9,:] = [0,+4]
		direct1[10,:] = [+4,0]
		direct1[11,:] = [0,-4]
		direct1[12,:] = [+2,-5]
		direct1[13,:] = [+2,+5]
		direct1[14,:] = [-2,+5]
		direct1[15,:] = [-2,-5]
		direct1[16,:] = [+5,-2]
		direct1[17,:] = [+5,+2]
		direct1[18,:] = [-5,+2]
		direct1[19,:] = [-5,-2]
		direct1[20,:] = [+2,-3]
		direct1[21,:] = [+2,+3]
		direct1[22,:] = [-2,+3]
		direct1[23,:] = [-2,-3]
		direct1[24,:] = [+3,-2]
		direct1[25,:] = [+3,+2]
		direct1[26,:] = [-3,+2]
		direct1[27,:] = [-3,-2]
		for i in range(num_mol):
			temp_count = corr_num(mol[i,:],direct1,latt)
			#if index[i] == 182:
			#	print "temp_count 108 : " + str(temp_count)
			#print "entry " + str(i) + "  temp_count " + str(temp_count)
			if temp_count > th:
				mol_count = mol_count + 1
				mol_list.append(index[i])
	return (mol_count,mol_list)
				
def cluster(mode, lattice):
	latt = deepcopy(lattice)
	latt = prepro(latt)
	mol = np.transpose(np.array(np.where(lattice == 3)))
	#print mol
	num_mol = mol.shape[0]
	output = np.array(range(1,num_mol+1))
	output_cluster = []
	if output.shape[0] != num_mol:
		print "WARNING!!! WRONG INITIALIZATION OF ELEMENT!!!"
	if mode == 0:
		direct = np.zeros((8,2))
		direct[0,:] = [-2,+1]
		direct[1,:] = [-1,+2]
		direct[2,:] = [+1,+2]
		direct[3,:] = [+2,+1]
		direct[4,:] = [+2,-1]
		direct[5,:] = [+1,-2]
		direct[6,:] = [-1,-2]
		direct[7,:] = [-2,-1]
	elif mode == 1:
		direct = np.zeros((12,2))
		direct[0,:] = [-2,+1]
		direct[1,:] = [-1,+2]
		direct[2,:] = [+1,+2]
		direct[3,:] = [+2,+1]
		direct[4,:] = [+2,-1]
		direct[5,:] = [+1,-2]
		direct[6,:] = [-1,-2]
		direct[7,:] = [-2,-1]
		direct[8,:] = [-2,0]
		direct[9,:] = [0,+2]
		direct[10,:] = [+2,0]
		direct[11,:] = [0,-2]
	elif mode == 2:
		direct = np.zeros((4,2))
		direct[0,:] = [-2,0]
		direct[1,:] = [0,+2]
		direct[2,:] = [+2,0]
		direct[3,:] = [0,-2]
	for i in range(0,num_mol):
		if output[i] == i+1:
			temp = [i] # list that store the connected elements
			templist = [i] # Queue, first in first out
			while len(templist) != 0:
				# pop out the first element
				current_index = templist.pop(0)
				# find out the neighbours of this element
				connected_neighbours = find_neighbour(mode,mol[current_index,:],direct,latt)
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
						#latt[mol[connected_neighbours[j]][0]][mol[connected_neighbours[j]][1]] = i+1			
			output_cluster.append(temp)
	print output_cluster
	new_mol = []
	num_count = 0
	for i in range(len(output_cluster)):
		if len(output_cluster[i]) > 3:
			ind = np.transpose(np.array(output_cluster[i]))
			ind_num = ind.shape[0]
			ele = mol[ind,:]
			count, mols = auto_correlate(mode, 2, ele, ind, lattice)
			if count > 3:
				new_mol.append(mols)
				num_count = num_count + len(mols)
	return new_mol, num_count
	#return (output_cluster,latt)

def get_coor_mol(input_coor):
	x = input_coor[0]
	y = input_coor[1]
	coor = np.zeros((5,2))
	coor[0,:] = [x, y]
	coor[1,:] = [sir(x-1,latt_len), y]
	coor[2,:] = [x, sir(y+1,latt_len)]
	coor[3,:] = [sir(x+1,latt_len), y]
	coor[4,:] = [x, sir(y-1,latt_len)]
	return coor

def cal_bond_num(latt):
	direct = np.zeros((4,2))
	direct[0,:] = [-1,0]
	direct[1,:] = [0,+1]
	direct[2,:] = [+1,0]
	direct[3,:] = [0,-1]
	cbond = np.zeros(5)
	temp1,temp2 = np.where(latt == 3)
	mol = np.array(zip(temp1,temp2))
	#print mol
	num_mol = mol.shape[0]
	#print num_mol
	temp1,temp2 = np.where(latt == 1)
	metal = np.array(zip(temp1,temp2))
	#Get the points around this molecule
	for i in range(0,num_mol):
		pos_around0 = get_coor_mol(mol[i])[1:5]
		pos_around1 = pos_around0 + direct
		pos_around2 = pos_around1 + direct
		pos_around3 = pos_around2 + direct
		count = 0
		for i in range(0,4):
			if latt[sir(pos_around1[i,0],latt_len),sir(pos_around1[i,1],latt_len)] == 1:
				if latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)] == 2 and \
						latt[sir(pos_around3[i,0],latt_len),sir(pos_around3[i,1],latt_len)] == 3:
							count = count + 1
			#elif latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != 0:
				#if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != \
						#latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)]: 
				#energy = energy - venergy
		cbond[count] = cbond[count] + 1
	return cbond

if __name__ == "__main__":
	dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results14"
	os.chdir(dname)
	latt_len = 80
	files = os.listdir(dname)
	line = files[32]
	print line
	lattice = np.loadtxt(line, delimiter=',',skiprows=1)
	mol = np.transpose(np.array(np.where(lattice == 3)))
	#print mol
	num_mol = mol.shape[0]
	cbond = cal_bond_num(lattice)
	print cbond
	output, count = cluster(1,lattice)
	for i in range(len(output)):
		if len(output[i]) > 3 :
			print "new: " + str(len(output[i])) + "  " + str(output[i])
	fig = plt.figure()	
	a = fig.add_subplot(1,2,1)
	imgplot = plt.imshow(lattice)
	new_latt1 = np.zeros((lattice.shape[0],lattice.shape[0]))
	print " count : " + str(count)
	for i in range(len(output)):
		if len(output[i]) > 3:
			for j in range(len(output[i])):
				#plt.text(mol[output[i][j]][1],mol[output[i][j]][0], str(output[i][j]),fontsize=8)
				new_latt1[mol[output[i][j]][0]][mol[output[i][j]][1]] = (i+1)*3
	a = fig.add_subplot(1,2,2)
	imgplot = plt.imshow(new_latt1)
	plt.show()
