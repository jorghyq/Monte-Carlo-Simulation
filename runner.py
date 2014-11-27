# runner class

import os
import subprocess
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt
import time


class McRunner:
	def __init__(self, mode):
		self.mode = mode
	
	def set_range_cenergy(self,cenergy_init,cenergy_step,cenergy_max):
		self.cenergy_init = cenergy_init
		self.cenergy_step = cenergy_step
		self.cenergy_max = cenergy_max		
		
	def set_range_nmetal(self,nmet_init,nmet_step,nmet_max):
		self.nmet_init = nmet_init
		self.nmet_step = nmet_step
		self.nmet_max = nmet_max
		
	def set_range_venergy(self,venergy_init,venergy_step,venergy_max):
		self.venergy_init = venergy_init
		self.venergy_step = venergy_step
		self.venergy_max = venergy_max
	
	def set_initial(self,total_run,num_mol,latt_len,folder_num):
		self.total_run = total_run
		self.num_mol = num_mol
		self.latt_len = latt_len
		self.folder_num = folder_num
		
	def run(self,EnableOutput):
		num_metal = self.nmet_init
		t1 = 0
		total_num = ((self.nmet_max - self.nmet_init)/self.nmet_step + 1) * \
			((self.cenergy_max - self.cenergy_init)/self.cenergy_step + 1) *\
			((self.venergy_max - self.venergy_init)/self.venergy_step + 1)
		print "total number of run :" + str(total_num)
		while num_metal <= self.nmet_max:
			cenergy = self.cenergy_init
			mcenergy = self.cenergy_init
			while cenergy <= self.cenergy_max:
				venergy = self.venergy_init
				while venergy <= self.venergy_max:
					os.system('mc-rect-lattice-func%d -a %d -b %d -c %d -d %d -e %f -f %d -g %d' % (self.mode,self.total_run,self.num_mol,num_metal,cenergy,venergy,mcenergy,self.folder_num))
					venergy = venergy + self.venergy_step
				cenergy = cenergy + self.cenergy_step
				mcenergy = mcenergy + self.cenergy_step	
			num_metal = num_metal + self.nmet_step
			if EnableOutput:
				t2 = time.clock()
				dt = t2 - t1
				t1 = t2
				print "This run costs time: %f" % (dt)
		
	def logfile(self,commet):
		f = open('D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results'+str(self.folder_num)+'\logfile.txt', 'w')
		f.write('This is the logfile for the following settings.\n')
		f.write(commet + '.\n')
		f.write('total_run: ' + str(self.total_run) + '\n')
		f.write('latt_len: ' + str(self.latt_len) + '\n')
		f.write('num_mol: ' + str(self.num_mol) + '\n')
		f.write('nmet_int: ' + str(self.nmet_init) + '\n')
		f.write('nmet_setp: ' + str(self.nmet_step) + '\n')
		f.write('nmet_max: ' + str(self.nmet_max) + '\n')
		f.write('cenergy_init: ' + str(self.cenergy_init) + '\n')
		f.write('cenergy_step: ' + str(self.cenergy_step) + '\n')
		f.write('cenergy_max: ' + str(self.cenergy_max) + '\n')
		f.write('venergy_init: ' + str(self.venergy_init) + '\n')		
		f.write('venergy_step: ' + str(self.venergy_step) + '\n')
		f.write('venergy_max: ' + str(self.venergy_max) + '\n')
		f.close()
		

if __name__ == "__main__":
	runner = McRunner(2)
	runner.set_initial(1000000,200,200,3)
	runner.set_range_cenergy(20,1,20)
	runner.set_range_venergy(0,2,30)
	runner.set_range_nmetal(50,50,450)
	runner.logfile('This is for the molecule BDS286')
	runner.run(1)	
