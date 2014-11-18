# Program to run the
import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time

total_run = 10000000
latt_len = 80
num_mol = 200
num_metal = 50
cenergy_initial = 30
venergy_initial = 3
mcenergy_initial = 30

num_metal_max = 410
cenergy_max = 60
num_metal_step = 50
cenergy_step = 2
###### generate the log file

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results5"
os.chdir(dname)
f = open('logfile-%d-%d-%d-%d-', 'w')
f.write('This is the logfile for the following settings.\n')
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
f.close()

while num_metal < num_metal_max:
	cenergy = cenergy_initial
	mcenergy = mcenergy_initial
	while cenergy < cenergy_max:
		os.system('mc-rect-lattice-func -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
		#lattice = np.loadtxt("results3\%.1e-%d-%d-%d-%d-%d.txt" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy), delimiter=',')
		#plt.imsave("results3\%.1e-%d-%d-%d-%d-%d.png" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy),lattice,[0,2])
		cenergy = cenergy + cenergy_step
		mcenergy = mcenergy + cenergy_step
		#print "num_metal = %d, venerg = %d" % (num_metal,venergy)
		
	num_metal = num_metal + num_metal_step
	

