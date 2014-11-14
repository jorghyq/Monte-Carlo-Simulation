# Make images from the text files
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt

def sir(x,range):
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
	coor[1,:] = [sir(x-1,latt_len), y]
	coor[2,:] = [x, sir(y+1,latt_len)]
	coor[3,:] = [sir(x+1,latt_len), y]
	coor[4,:] = [x, sir(y-1,latt_len)]
	return coor

def cal_bond_num(latt):
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
direct = np.zeros((4,2))
direct[0,:] = [-1,0]
direct[1,:] = [0,+1]
direct[2,:] = [+1,0]
direct[3,:] = [0,-1]

latt_len = 15

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results5"
os.chdir(dname)
files = os.listdir(dname)

xvalue = range(1,20,2)
xvalue2 = range(50,410,50)

cbond0 = np.zeros((10,8))
cbond1 = np.zeros((10,8))
cbond2 = np.zeros((10,8))
cbond3 = np.zeros((10,8))
cbond4 = np.zeros((10,8))

for line in files:
	if (line[-4:] == '.txt'):
		namedata = line[0:-4].strip().split('-')
		latt_len = int(namedata[1])
		num_metal = int(namedata[3])
		num_metal_ind = num_metal/50 - 1
		cenergy = int(namedata[4])
		cenergy_ind = (cenergy -1)/2
		#f = open(line, 'r')
		#headdata = f.readline().strip().split(',')
		#f.close()
		lattice = np.loadtxt(line, delimiter=',',skiprows=1)
		latt = lattice[0:latt_len,0:]
		bond_num = cal_bond_num(latt)
		cbond0[cenergy_ind][num_metal_ind] = bond_num[0]
		cbond1[cenergy_ind][num_metal_ind] = bond_num[1]
		cbond2[cenergy_ind][num_metal_ind] = bond_num[2]
		cbond3[cenergy_ind][num_metal_ind] = bond_num[3]
		cbond4[cenergy_ind][num_metal_ind] = bond_num[4]
		#lattice_num = lattice[latt_len:2*latt_len,0:]
		#cbondvalue[cenergy_ind][num_metal_ind] = float(headdata[0])/float(num_metal)
		#vbondvalue[cenergy_ind][num_metal_ind] = int(headdata[1])
		#energyvalue[cenergy_ind][num_metal_ind] = int(headdata[2])
print cbond0
cbond12 = cbond1 + cbond2
cbond34 = cbond3 + cbond4		
cbond4tr = np.transpose(cbond4)
print cbond4tr.shape
plt.figure()
for i in range(0,cbond34.shape[0]):
	plt.plot(xvalue2,cbond34[i,:],label = 'cenerg = %d' % (i*2+1))
plt.legend(loc=0,fontsize = 8)
plt.xlabel('metal numbers')
plt.ylabel('coordination bond numbers')
plt.show()
