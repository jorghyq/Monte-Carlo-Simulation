# call c++ many times with post process
# Program to run the
import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time
from analyzer import McAnalayzer

run_times = 10

total_run = 10000000
latt_len = 80
num_mol = 200
num_metal = 50
cenergy_initial = 1
venergy_initial = 1
mcenergy_initial = 1

num_metal_max = 450
cenergy_max = 29
num_metal_step = 50
cenergy_step = 2
venergy_max = 3
venergy_step = 2

nmetal_ind = (num_metal_max - num_metal)/num_metal_step + 1
cenergy_ind = (cenergy_max - num_metal)/cenergy_step + 1
venergy_ind = (venergy_max - venergy_initial)/venergy_step + 1

########################################################################
# 0.dim: venergy
# 1.dim: runtimes
# 2.dim: Ec
# 3.dim: nmetal
# 4.dim: total_energy, av_energy, total_coordination, av_coordination, dense, 1d, 2d, disordered
##########################################################################
results = np.zeros((run_times,cenergy_ind,nmetal_ind,8))

###### generate the log file

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results6"
os.chdir(dname)
f = open('D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results6\logfile.txt', 'w')
f.write('This is the logfile for the following settings.\n')
f.write('There is no constrained, no is_forbidden function.\n')
f.write('run_times: ' + str(run_times) + '\n')
f.write('total_run: ' + str(total_run) + '\n')
f.write('latt_len: ' + str(latt_len) + '\n')
f.write('num_mol: ' + str(num_mol) + '\n')
f.write('num_metal: ' + str(num_metal) + '\n')
f.write('cenergy_initial: ' + str(cenergy_initial) + '\n')
f.write('venergy_initial: ' + str(venergy_initial) + '\n')
f.write('mcenergy_initial: ' + str(mcenergy_initial) + '\n')
f.write('num_metal_max: ' + str(num_metal_max) + '\n')
f.write('cenergy_max: ' + str(cenergy_max) + '\n')
f.write('num_metal_step: ' + str(num_metal_step) + '\n')
f.write('cenergy_step: ' + str(cenergy_step) + '\n')
f.write('venergy_initial: ' + str(venergy_initial) + '\n')
f.write('venergy_step: ' + str(venergy_step) + '\n')
f.close()

az = McAnalayzer(dname)
analyzer.set_initial(nmetal_initial,nmetal_step,cenergy_initial,cenergy_step,venergy_initial,venergy_step)

#venergy = venergy_initial
while num_metal <= num_metal_max:
	cenergy = cenergy_initial
	mcenergy = mcenergy_initial
	while cenergy <= cenergy_max:
		venergy = venergy_initial
		while venergy <= venergy_max:
			for i in range(run_times):
				os.system('mc-rect-lattice-func2 -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
				
				str_venergy = '%.1e'% venergy + str(0)
				str_total_run = '1.0e+007'
				line = str_total_run+'-'+str(latt_len)+'-'+str(num_mol)+'-'+str(num_metal)+'-'+str(cenergy)+'-'+str_venergy+'-'+str(mcenergy)
				az.load_txt(line)
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][0] = az.totalenergy
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][1] = az.totalenergy_av
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][2] = az.cbond_num
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][3] = az.cbond_num_avt
				mols0, count0 = az.clustering(0)
				mols1, count1 = az.clustering(1)
				mols2, count2 = az.clustering(2)
				newlist = []
				count3 = 0
				mol_list = range(200)
				for i in range(len(mols0)):
					for j in range(len(mols0[i])):
						if mols0[i][j] not in newlist:
							newlist.append(mols0[i][j])
				for i in range(len(mols1)):
					for j in range(len(mols1[i])):
						if mols1[i][j] not in newlist:
							newlist.append(mols1[i][j])
				for i in range(len(mols2)):
					for j in range(len(mols2[i])):
						if mols2[i][j] not in newlist:
							newlist.append(mols2[i][j])
				for i in range(200):
					if i not in newlist:
						count3 = count3 + 1
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][4] = count0
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][5] = count1
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][6] = count2
				results[az.venergy_ind][i][az.cenergy_ind][az.nmetal_ind][7] = count3
				
				venergy = venergy + venergy_step
		#lattice = np.loadtxt("results3\%.1e-%d-%d-%d-%d-%d.txt" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy), delimiter=',')
		cenergy = cenergy + cenergy_step
		mcenergy = mcenergy + cenergy_step
		#print "num_metal = %d, venerg = %d" % (num_metal,venergy)
		
	num_metal = num_metal + num_metal_step
	

