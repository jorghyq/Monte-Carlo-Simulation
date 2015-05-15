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
		self.cenergy_init = float(cenergy_init)
		self.cenergy_step = float(cenergy_step)
		self.cenergy_max = float(cenergy_max)		
		
	def set_range_nmetal(self,nmet_init,nmet_step,nmet_max):
		self.nmet_init = nmet_init
		self.nmet_step = nmet_step
		self.nmet_max = nmet_max
		
	def set_range_venergy(self,venergy_init,venergy_step,venergy_max):
		self.venergy_init = float(venergy_init)
		self.venergy_step = float(venergy_step)
		self.venergy_max = float(venergy_max)
	
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
					os.system('./mc-rect-lattice-func-linux%d -a %d -b %d -c %d -d %f -e %f -f %f -g %d' % (self.mode,self.total_run,self.num_mol,num_metal,cenergy,venergy,mcenergy,self.folder_num))
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
		f.write('nmet_init: ' + str(self.nmet_init) + '\n')
		f.write('nmet_step: ' + str(self.nmet_step) + '\n')
		f.write('nmet_max: ' + str(self.nmet_max) + '\n')
		f.write('cenergy_init: ' + str(self.cenergy_init) + '\n')
		f.write('cenergy_step: ' + str(self.cenergy_step) + '\n')
		f.write('cenergy_max: ' + str(self.cenergy_max) + '\n')
		f.write('venergy_init: ' + str(self.venergy_init) + '\n')		
		f.write('venergy_step: ' + str(self.venergy_step) + '\n')
		f.write('venergy_max: ' + str(self.venergy_max) + '\n')
		f.close()
		

if __name__ == "__main__":
	runner = McRunner(1)
	runner.set_initial(10000000,300,100,37)
	#runner.set_range_cenergy(30,1,30)
	#runner.set_range_venergy(3,2,31)	runner.set_range_nmetal(0,25,300)
	#runner.logfile('This is for the molecule BDS285')
	#runner.run(1)
	#cenergy = 40
	#venergy = 40
	#kT = 3
	#while kT<14:
	#	app_cenergy = float(cenergy)/float(kT)
	#	app_venergy = float(app_cenergy/4)
	#	#app_venergy_step = float(app_venergy)/8
	#	runner.set_range_cenergy(app_cenergy,1,app_cenergy)
	#	runner.set_range_venergy(app_venergy,1,app_venergy)
	#	print "app_cenergy = " + str(app_cenergy) + " app_venergy = " + str(app_venergy)# +" app_venergy_step = " + str(app_venergy_step)
	#	runner.set_range_nmetal(0,50,350)
	#	runner.logfile('This is for the molecule BDS286')
	#	runner.run(1)
	#	kT = kT + 2
	kT_inv = 0.025
	while kT_inv <=0.3:
		app_cenergy = float(1/kT_inv)
		print "app_cenergy" + str(app_cenergy)
		venergy_step = float(app_cenergy/8)
		print "app_venergy_step" + str(venergy_step)
		runner.set_range_cenergy(app_cenergy,1,app_cenergy)
		runner.set_range_venergy(0,venergy_step,app_cenergy)
		runner.set_range_nmetal(300,1,300)
		runner.logfile('This is for the molecules BDS285')
		runner.run(1)
		kT_inv = kT_inv + 0.05	
