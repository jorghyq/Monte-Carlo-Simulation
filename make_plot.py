# make plot from data
import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results5"
os.chdir(dname)
latt_len = 80
files = os.listdir(dname)
xvalue = range(1,31,2)
xvalue2 = range(50,410,50)
#print xvalue
cbondvalue = np.zeros((15,8))
vbondvalue = np.zeros((15,8))
energyvalue = np.zeros((15,8))
for line in files:
	if (line[-4:] == '.txt'):
		namedata = line[0:-4].strip().split('-')
		num_metal = int(namedata[3])
		num_metal_ind = num_metal/50 - 1
		cenergy = int(namedata[4])
		cenergy_ind = (cenergy -1)/2
		
		#print namedata
		#print cenergy
		f = open(line, 'r')
		headdata = f.readline().strip().split(',')
		#latt_len = int(headdata[-1])
		#print headdata
		f.close()
		cbondvalue[cenergy_ind][num_metal_ind] = float(headdata[0])/float(num_metal)
		vbondvalue[cenergy_ind][num_metal_ind] = int(headdata[1])
		energyvalue[cenergy_ind][num_metal_ind] = int(headdata[2])

plt.figure()
#plt.plot(xvalue,cbondvalue)
#plt.plot(xvalue,vbondvalue)
#plt.plot(xvalue2,np.transpose(vbondvalue))
plt.plot(xvalue2,np.transpose(cbondvalue))
#plt.plot(xvalue,energyvalue)
plt.show()
