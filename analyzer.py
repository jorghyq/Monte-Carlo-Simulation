# Python file to analyze the lattice obtained from the simulation with the form specified below
# the center of the molecules : 3
# legs : 2
# metals : 1

import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster



def sprint(txt, inp):
	print txt + " : " + str(inp)


class McAnalyzer:
	def __init__(self, wd_path):
		self.path = wd_path
		
	def load_logfile(self):
		fname = os.path.join(self.path,'logfile.txt')
		if os.path.isfile(fname):
			f = open(fname, 'r')
			for line in f:
				tline = line.split('')
				if len(tline2) > 1:
					if tline2[0].strip() == 'total_run':
						self.total_run = float(tline2[1].strip())
					elif tline2[0].strip() == 'latt_len':
						self.latt_len = int(tline2[1].strip())	
					elif tline2[0].strip() == 'num_mol':
						self.num_mol = int(tline2[1].strip())
					elif tline2[0].strip() == 'nmet_init':
						self.nmet_init = float(tline2[1].strip())
					elif tline2[0].strip() == 'nmet_step':
						self.nmetal_step = float(tline2[1].strip())
					elif tline2[0].strip() == 'nmet_max':
						self.nmet_max = float(tline2[1].strip())
					elif tline2[0].strip() == 'cenergy_init':
						self.cenergy_init = float(tline2[1].strip())
					elif tline2[0].strip() == 'cenergy_step':
						self.cenergy_step = float(tline2[1].strip())
					elif tline2[0].strip() == 'cenergy_max':
						self.cenergy_max = float(tline2[1].strip())
					elif tline2[0].strip() == 'venergy_init':
						self.venergy_init = float(tline2[1].strip())
					elif tline2[0].strip() == 'venergy_step':
						self.venergy_step = float(tline2[1].strip())
					elif tline2[0].strip() == 'venergy_max':
						self.venergy_max = float(tline2[1].strip())
		else:
			print "no logfiles!!!!!"
		self.ind_metal = (self.nmet_max - nmet_init)/nmet_step + 1
		self.ind_cenergy = (self.cenergy_max - self.cenergy_init)/self.cenergy_step + 1
		self.ind_venergy = (self.venergy_max - self.venergy_init)/self.venergy_step + 1
		
	def set_initial(self,nmet_init,nmet_step,cenergy_init,cenergy_step,venergy_init,venergy_step):
		self.nmet_init = nmet_init
		self.nmet_step = nmet_step
		self.cenergy_init = cenergy_init
		self.cenergy_step = cenergy_step
		self.venergy_init = venergy_init
		self.venergy_step = venergy_step
		
	
	def load_txt(self,txt_name):
		#temp,txt_name = os.path.split(txt_name)
		####### information from the name #######
		temp, filename = os.path.split(txt_name)
		namedata = filename[0:-4].strip().split('-')
		self.latt_len = int(namedata[1])
		self.num_metal = int(namedata[3])
		#self.num_mol = int(namedata[2])
		self.nmetal_ind = (self.num_metal - self.nmet_init)/self.nmet_step
		self.cenergy = float(namedata[4])
		self.venergy = float(namedata[5])
		self.cenergy_ind = (self.cenergy - self.cenergy_init)/self.cenergy_step
		self.venergy_ind = (self.venergy - self.venergy_init)/self.venergy_step
		####### information from the header #######
		f = open(txt_name, 'r')
		headdata = f.readline().strip().split(',')
		f.close()
		self.cbond_num = float(headdata[0])
		self.cbond_num_av = self.cbond_num/float(self.num_metal)
		self.cbond_num_avt = self.cbond_num/float(self.num_metal + self.num_mol)	
		self.vbond_num = float(headdata[1])
		self.total_energy = float(headdata[2])
		self.energy_av = self.total_energy/float(self.num_metal + self.num_mol)
		####### information from the matrix #######
		self.lattice = np.loadtxt(txt_name, delimiter=',',skiprows=1)
		self.lattice = self.lattice[0:self.latt_len,0:]
		
	def run(self,mode):
		files = os.listdir(self.path)
		os.chdir(dname)
		if self.cenergy_int == self.cenergy_max:
			totalenergy = np.zeros((ind_metal,ind_venergy))
			cbond_num = np.zeros((ind_metal,ind_venergy))
			vbond_num = np.zeros((ind_metal,ind_venergy))
			cbond_num_av = np.zeros((ind_metal,ind_venergy))
			totalenergy_av = np.zeros((ind_metal,ind_venergy))
			cbond_num_avt = np.zeros((ind_metal,ind_venergy))
			mdense = np.zeros((ind_metal,ind_venergy))
			m1d = np.zeros((ind_metal,ind_venergy))
			m2d = np.zeros((ind_metal,ind_venergy))
			mdis = np.zeros((ind_metal,ind_venergy))
			for line in files:
				if line[0] == '1' and line[-4:] == '.txt':
					self.load_txt(line):
					totalenergy[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.total_energy
					cbond_num[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num
					vbond_num[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.vbond_num
					cbond_num_av[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num_av
					totalenergy_av[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.energy_av
					cbond_num_avt[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num_avt	
					if mode == 1:
						mols0, count0 = self.clustering(0)
						mdense[self.nmetal_ind][self.venergy_ind] = count0
						mols1, count1 = self.clustering(1)
						m1d[analyzer.nmetal_ind][self.venergy_ind] = count1
						mols2, count2 = self.clustering(2)
						m2d[self.nmetal_ind][self.venergy_ind] = count2
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
						mdis[self.nmetal_ind][self.venergy_ind] = count3
			np.savetxt("totalenergy.txt",totalenergy,delimiter=',')	
			np.savetxt("cbond_num.txt",cbond_num,delimiter=',')
			np.savetxt("vbond_num.txt",vbond_num,delimiter=',')
			np.savetxt("cbond_num_av.txt",cbond_num_av,delimiter=',')
			np.savetxt("totalenergy_av.txt",totalenergy_av,delimiter=',')
			np.savetxt("cbond_num_avt.txt",cbond_num_avt,delimiter=',')
			np.savetxt("mdense.txt",mdense, delimiter=',')	
			np.savetxt("m1d.txt",m1d, delimiter=',')
			np.savetxt("m2d.txt",m2d, delimiter=',')
			np.savetxt("mdis.txt",mdis, delimiter=',')
	
	def plot(filename, mode,xlab, ylab):
		fname = os.path.join(self.path,filename)
		temp_file = np.loadtxt(fname,delimiter=',')
		if os.path.exists(temp_file):
			if mode == 1:
				x_axis = np.array(range(cenergy_init,venergy_max+venergy_step,venergy_step))
				for i in range(0,temp_file.shape[1]):
					plt.plot(x_axis,temp_file[:,i],label = 'cenerg = %d' % (i*2+1))
					plt.legend(loc=0,fontsize = 8)
					plt.xlabel(xlab)
					plt.ylabel(ylab)
		else:
			print "no input file!"
	
	def bond_num(self):
		temp = cal_bond_num(self.lattice) # 1*5 vector
		return temp
	
	def clustering(self,mode):
		mols,count = cluster(mode,self.lattice) # mols: list, count: total number
		return mols,count 

#class Summarizer:
#	def __init__(self):

if __name__ == "__main__":
	# go to the working directory
	dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results2"

	#######################################
	# dmol_num: ind_metal * ind_cenergy
	#   1) dense
	#   2) 1d
	#   3) 2d
	#	4) disordered
	########################################	
	np.savetxt("totalenergy.txt",totalenergy,delimiter=',')	
	np.savetxt("cbond_num.txt",cbond_num,delimiter=',')
	np.savetxt("vbond_num.txt",vbond_num,delimiter=',')
	np.savetxt("cbond_num_av.txt",cbond_num_av,delimiter=',')
	np.savetxt("totalenergy_av.txt",totalenergy_av,delimiter=',')
	np.savetxt("cbond_num_avt.txt",cbond_num_avt,delimiter=',')
	
	#np.savetxt("mdense.txt",mdense, delimiter=',')	
	#np.savetxt("m1d.txt",m1d, delimiter=',')
	#np.savetxt("m2d.txt",m2d, delimiter=',')
	#np.savetxt("mdis.txt",mdis, delimiter=',')
	#mdense = np.loadtxt("mdense.txt", delimiter=',')
	#m1d = np.loadtxt("m1d.txt", delimiter=',')
	#m2d = np.loadtxt("m2d.txt", delimiter=',')
	#mdis = np.loadtxt("mdis.txt", delimiter=',')
	#mtotal = mdense + m1d + m2d + mdis
	#totalenergy_av = np.loadtxt("totalenergy_av.txt",delimiter=',')
	#cbond_num_avt = np.loadtxt("cbond_num_avt.txt",delimiter=',')
	
	
	#mdense_p = mdense/mtotal
	#m1d_p = m1d/mtotal
	#m2d_p = m2d/mtotal
	#mdis_p = mdis/mtotal
	xmetal = np.array(range(nmetal_initial,nmetal_final + 10,nmetal_step))
	##print xmetal.shape
	#fig = plt.figure()
	#fig.add_subplot(2,2,1)
	#for i in range(0,mdense.shape[1]):
		#plt.plot(xmetal,mdense_p[:,i],label = 'cenerg = %d' % (i*2+1))
		#plt.text(xmetal[int(ind_metal/2)],mdense_p[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	#plt.legend(loc=0,fontsize = 8)
	#plt.xlabel('number of metals')
	#plt.ylabel('number of molecules in dense-packed')
	
	
	#fig.add_subplot(2,2,2)
	#for i in range(0,m1d.shape[1]):
		#plt.plot(xmetal,m1d_p[:,i],label = 'cenerg = %d' % (i*2+1))
		#plt.text(xmetal[int(ind_metal/2)],m1d_p[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	#plt.legend(loc=0,fontsize = 8)
	#plt.xlabel('number of metals')
	#plt.ylabel('number of molecules in 1d')
	
	#fig.add_subplot(2,2,3)
	#for i in range(0,m2d.shape[1]):
		#plt.plot(xmetal,m2d_p[:,i],label = 'cenerg = %d' % (i*2+1))
		#plt.text(xmetal[int(ind_metal/2)],m2d_p[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	#plt.legend(loc=0,fontsize = 8)
	#plt.xlabel('number of metals')
	#plt.ylabel('number of molecules in 2d')
	
	#fig.add_subplot(2,2,4)
	#for i in range(0,mdis.shape[1]):
		#plt.plot(xmetal,mdis_p[:,i],label = 'cenerg = %d' % (i*2+1))
		#plt.text(xmetal[int(ind_metal/2)],mdis_p[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	#plt.legend(loc=0,fontsize = 8)
	#plt.xlabel('number of metals')
	#plt.ylabel('number of molecules in disordered')
	#plt.show()
	
	
	
	##################################################################
	fig = plt.figure()
	fig.add_subplot(2,3,1)
	for i in range(0,totalenergy.shape[1]):
		plt.plot(xmetal,totalenergy[:,i],label = 'cenerg = %d' % (i*2+1))
		plt.text(xmetal[int(ind_metal/2)],totalenergy[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	plt.legend(loc=0)
	plt.xlabel('number of metals')
	plt.ylabel('total energy')
	
	
	fig.add_subplot(2,3,2)
	for i in range(0,cbond_num.shape[1]):
		plt.plot(xmetal,cbond_num[:,i],label = 'cenerg = %d' % (i*2+1))
		plt.text(xmetal[int(ind_metal/2)],cbond_num[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	plt.legend(loc=0)
	plt.xlabel('number of metals')
	plt.ylabel('number of coordination bond')
	
	fig.add_subplot(2,3,3)
	for i in range(0,vbond_num.shape[1]):
		plt.plot(xmetal,vbond_num[:,i],label = 'cenerg = %d' % (i*2+1))
		plt.text(xmetal[int(ind_metal/2)],vbond_num[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	plt.legend(loc=0)
	plt.xlabel('number of metals')
	plt.ylabel('number of vdW bond')
	
	fig.add_subplot(2,3,4)
	for i in range(0,cbond_num_av.shape[1]):
		plt.plot(xmetal,cbond_num_av[:,i],label = 'cenerg = %d' % (i*2+1))
		plt.text(xmetal[int(ind_metal/2)],cbond_num_av[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	plt.legend(loc=0)
	plt.xlabel('number of metals')
	plt.ylabel('total coordination bond per metal')
	
	fig.add_subplot(2,3,5)
	for i in range(0,totalenergy_av.shape[1]):
		plt.plot(xmetal,totalenergy_av[:,i],label = 'cenerg = %d' % (i*2+1))
		plt.text(xmetal[int(ind_metal/2)],totalenergy_av[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	plt.legend(loc=0)
	plt.xlabel('number of metals')
	plt.ylabel('total energy per element')
	
	fig.add_subplot(2,3,6)
	for i in range(0,cbond_num_avt.shape[1]):
		plt.plot(xmetal,cbond_num_avt[:,i],label = 'cenerg = %d' % (i*2+1))
		plt.text(xmetal[int(ind_metal/2)],cbond_num_avt[int(ind_metal/2),i], str(i*2+1),fontsize=8)
	plt.legend(loc=0)
	plt.xlabel('number of metals')
	plt.ylabel('total coordination bond per element')
	plt.show()
	
