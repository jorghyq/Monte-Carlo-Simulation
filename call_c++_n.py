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
from analyzer import McAnalyzer

run_times = 5

total_run = 10000000
latt_len = 100
num_mol = 200
num_metal = 50
cenergy_initial = 1
venergy_initial = 1
mcenergy_initial = 1

num_metal_max = 450
cenergy_max = 35
num_metal_step = 50
cenergy_step = 2
venergy_max = 3
venergy_step = 2

nmetal_ind = (num_metal_max - num_metal)/num_metal_step + 1
cenergy_ind = (cenergy_max - cenergy_initial)/cenergy_step + 1
venergy_ind = (venergy_max - venergy_initial)/venergy_step + 1
print nmetal_ind
print cenergy_ind
print venergy_ind
########################################################################
# 0.dim: venergy
# 1.dim: runtimes
# 2.dim: Ec
# 3.dim: nmetal
# 4.dim: total_energy, av_energy, total_coordination, av_coordination, dense, 1d, 2d, disordered
##########################################################################
results = np.zeros((venergy_ind,run_times,cenergy_ind,nmetal_ind,8))

###### generate the log file

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results3"
#os.chdir(dname)
f = open('D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results3\logfile.txt', 'w')
f.write('This is the logfile for the following settings.\n')
#f.write('There is no constrained, no is_forbidden function.\n')
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

az = McAnalyzer(dname)
az.set_initial(num_metal,num_metal_step,cenergy_initial,cenergy_step,venergy_initial,venergy_step)



#venergy = venergy_initial
while num_metal <= num_metal_max:
	cenergy = cenergy_initial
	mcenergy = mcenergy_initial
	while cenergy <= cenergy_max:
		venergy = venergy_initial
		while venergy <= venergy_max:
			i = 0
			while i < run_times:
				os.system('mc-rect-lattice-func -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
				str_venergy = '%.1e'% venergy + str(0)
				str_total_run = '1.0e+007'
				line = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results3\\"+str_total_run+'-'+str(latt_len)+'-'+str(num_mol)+'-'+str(num_metal)+'-'+str(cenergy)+'-'+str_venergy+'-'+str(mcenergy)+'.txt'
				print line
				#tpath,line = os.path.split(line)
				#print line
				az.load_txt(line)
				print int(az.cenergy_ind)
				print int(az.venergy_ind)
				print int(az.nmetal_ind)
				

				
				mols0, count0 = az.clustering(0)
				mols1, count1 = az.clustering(1)
				mols2, count2 = az.clustering(2)
				newlist = []
				count3 = 0
				mol_list = range(200)
				for j in range(len(mols0)):
					for k in range(len(mols0[j])):
						if mols0[j][k] not in newlist:
							newlist.append(mols0[j][k])
				for j in range(len(mols1)):
					for k in range(len(mols1[j])):
						if mols1[j][k] not in newlist:
							newlist.append(mols1[j][k])
				for j in range(len(mols2)):
					for k in range(len(mols2[j])):
						if mols2[j][k] not in newlist:
							newlist.append(mols2[j][k])
				for j in range(200):
					if j not in newlist:
						count3 = count3 + 1
				
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][0] = az.total_energy
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][1] = az.energy_av
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][2] = az.cbond_num
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][3] = az.cbond_num_avt
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][4] = count0
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][5] = count1
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][6] = count2
				results[int(az.venergy_ind)][i][int(az.cenergy_ind)][int(az.nmetal_ind)][7] = count3
				
				print "num_mtetal: " + str(num_metal) + " run: " + str(i)  + " cenergy : " + str(cenergy) + " venergy: " + str(venergy)
				i = i + 1
			venergy = venergy + venergy_step
		#lattice = np.loadtxt("results3\%.1e-%d-%d-%d-%d-%d.txt" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy), delimiter=',')
		cenergy = cenergy + cenergy_step
		mcenergy = mcenergy + cenergy_step
		#print "num_metal = %d, venerg = %d" % (num_metal,venergy)
		
	num_metal = num_metal + num_metal_step


#print results 

print "Done..."
np.save("D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results3\\results",results)
#np.savetxt("results.txt",results,delimiter=',',fmt='%1.3f')	
