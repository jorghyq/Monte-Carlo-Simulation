# Program to run the
import os
import sys
import subprocess
import shlex
#from subprocess import call
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time
import json

total_run = 1000000
ffn=100
latt_len = 100
restore = 1

#print total_run, num_mol1, num_mol2, num_metal
cenergy_init = 40
cenergy = cenergy_init
cenergy_max = 40
cenergy_step = 1

mcenergy_init = 40
mcenergy = mcenergy_init
mcenergy_max = 40
mcenergy_step = 1

venergy_init = 5
venergy = venergy_init
venergy_max = 5
venergy_step = 1

num_mol1_init = 150
num_mol1 = num_mol1_init
num_mol1_max = 150
num_mol1_step = 1

num_mol2_init = 50
num_mol2 = num_mol2_init
num_mol2_max = 200
num_mol2_step = 50

num_metal_init = 50
num_metal = num_metal_init
num_metal_max = 300
num_metal_step = 50

while venergy <= venergy_max:
    while cenergy <= cenergy_max:
        while mcenergy <= mcenergy_max:
            while num_mol1 <= num_mol1_max:
                while num_mol2 <= num_mol2_max:
                    while num_metal <= num_metal_max:
                        command = './mc-rect-lattice-func-linux6 -a %d -b %d -c %d -d %d -e %f -f %f -g %f -h %d -i %d' % (total_run,num_mol1,num_mol2,num_metal,cenergy,venergy,mcenergy,ffn,restore)
                        args = shlex.split(command)
                        #print args
                        os.system(command)
                        num_metal = num_metal + num_metal_step
                        #print venergy, cenergy, mcenergy, num_mol1, num_mol2, num_metal
                    num_metal = num_metal_init
                    num_mol2 = num_mol2 + num_mol2_step
                num_metal = num_metal_init
                num_mol2 = num_mol2_init
                num_mol1 = num_mol1 + num_mol1_step
            num_metal = num_metal_init
            num_mol2 = num_mol2_init
            num_mol1 = num_mol1_init
            mcenergy = mcenergy + mcenergy_step
        num_metal = num_metal_init
        num_mol2 = num_mol2_init
        num_mol1 = num_mol1_init
        mcenergy = mcenergy_init
        cenergy = cenergy + cenergy_step
    num_metal = num_metal_init
    num_mol2 = num_mol2_init
    num_mol1 = num_mol1_init
    mcenergy = mcenergy_init
    cenergy = cenergy_init
    venergy = venergy + venergy_step
