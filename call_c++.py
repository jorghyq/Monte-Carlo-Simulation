# Program to run the
import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time

# weired
total_run = 10000000
latt_len = 100
num_mol = 200
num_metal = 50
cenergy_initial = 5
venergy_initial = 10
mcenergy_initial = 5

num_metal_max = 450
cenergy_max = 90
num_metal_step = 50

venergy_step = 1
cenergy_step = 8
venergy_max = 10
###### generate the log file

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results9"
#os.chdir(dname)
f = open('D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results9\logfile.txt', 'w')
f.write('This is the logfile for the following settings.\n')
f.write('There is no constrained, no is_forbidden function.\n')
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

venergy = venergy_initial
while num_metal <= num_metal_max:
	cenergy = cenergy_initial
	mcenergy = mcenergy_initial
	while cenergy <= cenergy_max:
		venergy = venergy_initial
		while venergy <= venergy_max:
			os.system('mc-rect-lattice-func2 -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
			venergy = venergy + venergy_step
		#lattice = np.loadtxt("results3\%.1e-%d-%d-%d-%d-%d.txt" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy), delimiter=',')
		#plt.imsave("results3\%.1e-%d-%d-%d-%d-%d.png" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy),lattice,[0,2])
		cenergy = cenergy + cenergy_step
		mcenergy = mcenergy + cenergy_step
		#print "num_metal = %d, venerg = %d" % (num_metal,venergy)
		
	num_metal = num_metal + num_metal_step
	
#while num_metal <= num_metal_max:
	#venergy = venergy_initial
	#while venergy <= venergy_max:
		#cenergy = venergy
		#mcenergy = venergy
		#for i in range(10):
			#os.system('mc-rect-lattice-func2 -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
			#cenergy = cenergy + cenergy_step
			#mcenergy = mcenergy + cenergy_step
		#venergy = venergy + venergy_step
		##lattice = np.loadtxt("results3\%.1e-%d-%d-%d-%d-%d.txt" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy), delimiter=',')
		##plt.imsave("results3\%.1e-%d-%d-%d-%d-%d.png" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy),lattice,[0,2])
		
		##print "num_metal = %d, venerg = %d" % (num_metal,venergy)
	#num_metal = num_metal + num_metal_step
