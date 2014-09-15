# Make images from the text files
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results6"
os.chdir(dname)
#print dname
#files  = os.listdir(dname)
#print files
#out = open("sta_date.txt","w")
files = os.listdir(dname)
for line in files:
	if (line[-4:] == '.txt'):
		newline = line[:-4]+'.png'
		if newline not in files:
			lattice = np.loadtxt(line, delimiter=',')
			#cbond, vbond = cal_bond_num(lattice)
			#out.write(line[-4:]+'\t')
			plt.savefig(newline)
			print newline
			plt.imsave(newline,lattice,[0,2])
	

#def cal_bond_num(latt):
    ##First detect if there are molecules around
    #cbond = 0
    #vbond = 0
    #temp1, temp2 = np.where(latt == 1)
    #temp = get_coor_mol(coor)
    #direct = np.zeros((4,2))
    #direct[0,:] = [-1,0]
    #direct[1,:] = [0,+1]
    #direct[2,:] = [+1,0]
    #direct[3,:] = [0,-1]
    ##Get the points around this molecule
    #pos_around = np.zeros((4,2))
    #pos_around = temp[1:5,:] + direct
    #pos_around2 = pos_around + direct
    #for i in range(0,4):
        #if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] == num_mol:
            #energy = energy - cenergy
        #elif latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != 0:
            #if latt[sir(pos_around[i,0],latt_len),sir(pos_around[i,1],latt_len)] != \
                    #latt[sir(pos_around2[i,0],latt_len),sir(pos_around2[i,1],latt_len)]: 
			#energy = energy - venergy
    #return energy
