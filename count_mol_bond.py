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
line = files[31]
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

def find_neighbour(mode, coor, direct, lattice):
	
	#print "input coordinate for find_neighbour: " + str(coor)
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
				
def cluster(mode, lattice, direct):
	
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
	# 1. cluster
	#output,latt = cluster(lattice)
	# 2. auto-correlation
	latt_len = lattice.shape[0]
	latt = np.zeros((latt_len,latt_len))
	#latt[tuple(map(tuple,element))] == 3
	for i in range(element.shape[0]):
		latt[tuple(element[i,:])] = 3
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


mol = np.transpose(np.array(np.where(lattice == 3)))
num_mol = mol.shape[0]
		
fig = plt.figure()	
a = fig.add_subplot(2,2,1)
imgplot = plt.imshow(lattice)



############################### Dense packed ############################
a = fig.add_subplot(2,2,2)

direct0 = np.zeros((8,2))
direct0[0,:] = [-2,+1]
direct0[1,:] = [-1,+2]
direct0[2,:] = [+1,+2]
direct0[3,:] = [+2,+1]
direct0[4,:] = [+2,-1]
direct0[5,:] = [+1,-2]
direct0[6,:] = [-1,-2]
direct0[7,:] = [-2,-1]

output0,latt0 = cluster(0, lattice, direct0)
new_mol0 = []

print "############# Here begins the dense packed information ###########"
for i in range(len(output0)):
	if len(output0[i]) > 3 :
		print "new: " + str(len(output0[i])) + "  " + str(output0[i])

#print "##################################################"
for i in range(len(output0)):
	if len(output0[i]) > 3:
		ind0 = np.transpose(np.array(output0[i]))
		ind_num0 = ind0.shape[0]
		ele0 = mol[ind0,:]
		count0, mols0 = auto_correlate(0, 2, ele0, ind0, lattice)
		if count0 > 3:
			new_mol0.append(mols0)

for i in range(len(new_mol0)):
	print "new: " + str(len(new_mol0[i])) + "  " + str(new_mol0[i])

new_latt0 = np.zeros((lattice.shape[0],lattice.shape[0]))


#imgplot = plt.imshow(lattice)
for i in range(len(new_mol0)):
	for j in range(len(new_mol0[i])):
		plt.text(mol[new_mol0[i][j]][1],mol[new_mol0[i][j]][0], str(new_mol0[i][j]),fontsize=10)
		new_latt0[mol[new_mol0[i][j]][0]][mol[new_mol0[i][j]][1]] = (i+1)*3

#for i in range(len(output0)):
#	if len(output0[i]) > 3:
#		for j in range(len(output0[i])):
#			plt.text(mol[output0[i][j]][1],mol[output0[i][j]][0], str(output0[i][j]),fontsize=8)
#			new_latt0[mol[output0[i][j]][0]][mol[output0[i][j]][1]] = (i+1)*3

imgplot = plt.imshow(new_latt0)

############################## 1D networks ################################

a = fig.add_subplot(2,2,3)

direct1 = np.zeros((12,2))
direct1[0,:] = [-2,+1]
direct1[1,:] = [-1,+2]
direct1[2,:] = [+1,+2]
direct1[3,:] = [+2,+1]
direct1[4,:] = [+2,-1]
direct1[5,:] = [+1,-2]
direct1[6,:] = [-1,-2]
direct1[7,:] = [-2,-1]
direct1[8,:] = [-2,0]
direct1[9,:] = [0,+2]
direct1[10,:] = [+2,0]
direct1[11,:] = [0,-2]

output1,latt1 = cluster(1, lattice, direct1)
new_mol1 = []
new_latt1 = np.zeros((lattice.shape[0],lattice.shape[0]))

print "############# Here begins the 1D information ###########"
for i in range(len(output1)):
	if len(output1[i]) > 3 :
		print "new: " + str(len(output1[i])) + "  " + str(output1[i])

#for i in range(len(output1)):
#	if len(output1[i]) > 3:
#		for j in range(len(output1[i])):
#			plt.text(mol[output1[i][j]][1],mol[output1[i][j]][0], str(output1[i][j]),fontsize=8)
#			new_latt1[mol[output1[i][j]][0]][mol[output1[i][j]][1]] = (i+1)*3
			
for i in range(len(output1)):
	if len(output1[i]) > 3:
		ind1 = np.transpose(np.array(output1[i]))
		ind_num1 = ind1.shape[0]
		ele1 = mol[ind1,:]
		count1, mols1 = auto_correlate(1, 2, ele1, ind1, lattice)
		if count1 > 3:
			new_mol1.append(mols1)
print "###############################"
for i in range(len(new_mol1)):
	print "new: " + str(len(new_mol1[i])) + "  " + str(new_mol1[i])
	
for i in range(len(new_mol1)):
	for j in range(len(new_mol1[i])):
		plt.text(mol[new_mol1[i][j]][1],mol[new_mol1[i][j]][0], str(new_mol1[i][j]),fontsize=10)
		#print new_mol2[i][j]
		new_latt1[mol[new_mol1[i][j]][0]][mol[new_mol1[i][j]][1]] = (i+1)*3		

imgplot = plt.imshow(new_latt1)






############################## 1D END ################################
################################ 2D networks ##############################

a = fig.add_subplot(2,2,4)
direct2 = np.zeros((4,2))
direct2[0,:] = [-2,0]
direct2[1,:] = [0,+2]
direct2[2,:] = [+2,0]
direct2[3,:] = [0,-2]
output2,latt2 = cluster(2,lattice, direct2)

print "############# Here begins 2D information ###########"
for i in range(len(output2)):
	if len(output2[i]) > 3 :
		print "new: " + str(len(output2[i])) + "  " + str(output2[i])


new_mol2 = []
for i in range(len(output2)):
	if len(output2[i]) > 3:
		ind2 = np.transpose(np.array(output2[i]))
		ind_num2 = ind2.shape[0]
		ele2 = mol[ind2,:]
		count2, mols2 = auto_correlate(2, 2, ele2, ind2, lattice)
		new_mol2.append(mols2)
print "###############################"
for i in range(len(new_mol2)):
	print "new: " + str(len(new_mol2[i])) + "  " + str(new_mol2[i])


new_latt2 = np.zeros((lattice.shape[0],lattice.shape[0]))
for i in range(len(new_mol2)):
	for j in range(len(new_mol2[i])):
		plt.text(mol[new_mol2[i][j]][1],mol[new_mol2[i][j]][0], str(new_mol2[i][j]),fontsize=10)
		#print new_mol2[i][j]
		new_latt2[mol[new_mol2[i][j]][0]][mol[new_mol2[i][j]][1]] = (i+1)*3	
		
imgplot = plt.imshow(new_latt2)
############################## 2D Networks END ###########################
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

