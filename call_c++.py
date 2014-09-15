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
cenergy = 10
venergy = 3
mcenergy = 10	
while num_metal < 410:
	cenergy = 3
	mcenergy = 3
	while cenergy < 30:
		os.system('mc-rect-lattice-func -a %d -b %d -c %d -d %d -e %f -f %d' % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy))
		#lattice = np.loadtxt("results3\%.1e-%d-%d-%d-%d-%d.txt" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy), delimiter=',')
		#plt.imsave("results3\%.1e-%d-%d-%d-%d-%d.png" % (total_run,num_mol,num_metal,cenergy,venergy,mcenergy),lattice,[0,2])
		cenergy = cenergy + 1
		mcenergy = mcenergy + 1
		#print "num_metal = %d, venerg = %d" % (num_metal,venergy)
		
	num_metal = num_metal + 50
