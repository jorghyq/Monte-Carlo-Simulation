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
from PIL import Image


def sprint(txt, inp):
	print txt + " : " + str(inp)


class McAnalyzer:
	def __init__(self, wd_path):
		self.path = wd_path
		
	def load_logfile(self):
		fname = os.path.join(self.path,'logfile.txt')
		if os.path.isfile(fname):
			f = open(fname, 'r')
			print "load logfile!"
			for line in f:
				tline = line.split(':')
				if len(tline) > 1:
					if tline[0].strip() == 'total_run':
						self.total_run = float(tline[1].strip())
					elif tline[0].strip() == 'latt_len':
						self.latt_len = int(tline[1].strip())	
					elif tline[0].strip() == 'num_mol':
						self.num_mol = int(tline[1].strip())
					elif tline[0].strip() == 'nmet_init':
						self.nmet_init = float(tline[1].strip())
					elif tline[0].strip() == 'nmet_step':
						self.nmet_step = float(tline[1].strip())
					elif tline[0].strip() == 'nmet_max':
						self.nmet_max = float(tline[1].strip())
					elif tline[0].strip() == 'cenergy_init':
						self.cenergy_init = float(tline[1].strip())
					elif tline[0].strip() == 'cenergy_step':
						self.cenergy_step = float(tline[1].strip())
					elif tline[0].strip() == 'cenergy_max':
						self.cenergy_max = float(tline[1].strip())
					elif tline[0].strip() == 'venergy_init':
						self.venergy_init = float(tline[1].strip())
					elif tline[0].strip() == 'venergy_step':
						self.venergy_step = float(tline[1].strip())
					elif tline[0].strip() == 'venergy_max':
						self.venergy_max = float(tline[1].strip())
		else:
			print "no logfiles!!!!!"
		print "logfile loaded!"
		self.ind_metal = (self.nmet_max - self.nmet_init)/self.nmet_step + 1
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
		#self.cbond_num_av = self.cbond_num/float(self.num_metal)
		self.cbond_num_avt = self.cbond_num/float(self.num_metal + self.num_mol)	
		self.vbond_num = float(headdata[1])
		self.total_energy = float(headdata[2])
		self.energy_av = self.total_energy/float(self.num_metal + self.num_mol)
		####### information from the matrix #######
		self.lattice = np.loadtxt(txt_name, delimiter=',',skiprows=1)
		self.lattice = self.lattice[0:self.latt_len,0:]
		
	def run(self,mode):
		files = os.listdir(self.path)
		#mdense = np.loadtxt("mdense.txt", delimiter=',')
		#m1d = np.loadtxt("m1d.txt", delimiter=',')
		#m2d = np.loadtxt("m2d.txt", delimiter=',')
		#mdis = np.loadtxt("mdis.txt", delimiter=',')
		#mtotal = mdense + m1d + m2d + mdis
		#fname = os.path.join(self.path,filename)
		os.chdir(dname)
		if self.cenergy_init == self.cenergy_max:
			totalenergy = np.zeros((self.ind_metal,self.ind_venergy))
			cbond_num = np.zeros((self.ind_metal,self.ind_venergy))
			vbond_num = np.zeros((self.ind_metal,self.ind_venergy))
			cbond_num_av = np.zeros((self.ind_metal,self.ind_venergy))
			totalenergy_av = np.zeros((self.ind_metal,self.ind_venergy))
			cbond_num_avt = np.zeros((self.ind_metal,self.ind_venergy))
			mdense = np.zeros((self.ind_metal,self.ind_venergy))
			m1d = np.zeros((self.ind_metal,self.ind_venergy))
			m2d = np.zeros((self.ind_metal,self.ind_venergy))
			mdis = np.zeros((self.ind_metal,self.ind_venergy))
			if os.path.exists(os.path.join(self.path,"mdense.txt")):
				mdense = np.loadtxt("mdense.txt", delimiter=',')
				if os.path.exists(os.path.join(self.path,"mdense.txt")):
					m1d = np.loadtxt("m1d.txt", delimiter=',')
					if os.path.exists(os.path.join(self.path,"mdense.txt")):
						m2d = np.loadtxt("m2d.txt", delimiter=',')
						if os.path.exists(os.path.join(self.path,"mdense.txt")):
							mdis = np.loadtxt("mdis.txt", delimiter=',')
							print "DATA ALREADY EXIST!"
							return
			print "Begin to process the data..."
			for line in files:
				
				if line[0] == '1' and line[-4:] == '.txt':
					self.load_txt(line)
					totalenergy[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.total_energy
					cbond_num[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num
					vbond_num[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.vbond_num
					#cbond_num_av[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num_av
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
			mtotal = mdense + m1d + m2d + mdis
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
			np.savetxt("mtotal.txt",mtotal, delimiter=',')
			print "Done, data are saved!"
	
	def plot_curve(self,filename, mode,xlab, ylab, prozent):
		fname = os.path.join(self.path,filename)
		temp_file = np.loadtxt(fname,delimiter=',')
		#print temp_file.shape
		temp_total = np.loadtxt('mtotal.txt',delimiter=',')
		if prozent == 1:
			temp_file = temp_file/temp_total
		if os.path.exists(fname):
			if mode == 1:
				x_axis = np.array(range(int(self.ind_metal)))
				#print int(self.ind_venergy)
				for i in range(0,temp_file.shape[1]):
					plt.plot(x_axis,temp_file[:,i],label = 'cenerg = %d' % (i*2+1))
					plt.legend(loc=0,prop={'size':8})
					plt.xlabel(xlab)
					plt.ylabel(ylab)
		else:
			print "no input file!"
	
	def phase_diagram(self,updown,leftright,xlab,ylab):
		mdense = np.loadtxt("mdense.txt", delimiter=',')
		m1d = np.loadtxt("m1d.txt", delimiter=',')
		m2d = np.loadtxt("m2d.txt", delimiter=',')
		mdis = np.loadtxt("mdis.txt", delimiter=',')
		mtotal = np.loadtxt('mtotal.txt',delimiter=',')
		mdense_p = mdense/mtotal
		m1d_p = m1d/mtotal
		m2d_p = m2d/mtotal
		if updown:
			mdense_p = np.flipud(mdense_p)
			m1d_p = np.flipud(m1d_p)
			m2d_p = np.flipud(m2d_p)
		if leftright:
			mdense_p = np.fliplr(mdense_p)
			m1d_p = np.fliplr(m1d_p)
			m2d_p = np.fliplr(m2d_p)
		r = m1d_p
		g = m2d_p
		b = mdense_p
		rgb = np.dstack((r,g,b))
		im = Image.fromarray(np.uint8(rgb*255.999))
		plt.imshow(im,extent=[0.1,2,50,450],aspect="auto")
		plt.xlabel(xlab)
		plt.ylabel(ylab)

	def bond_num(self):
		temp = cal_bond_num(self.lattice) # 1*5 vector
		return temp
	
	def clustering(self,mode):
		mols,count = cluster(mode,self.lattice) # mols: list, count: total number
		return mols,count 


if __name__ == "__main__":
	# go to the working directory
	dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results4"
	analyzer = McAnalyzer(dname)
	analyzer.load_logfile()
	analyzer.run(1)
	fig = plt.figure()
	analyzer.phase_diagram(0,0,"Ev/Ec","nmetal")
	#fig.add_subplot(2,2,1)
	#analyzer.plot_curve("mdense.txt",1,"number metals","prozent",1)

	#fig.add_subplot(2,2,2)
	#analyzer.plot_curve("m1d.txt",1,"number metals","prozent",1)
	
	#fig.add_subplot(2,2,3)
	#analyzer.plot_curve("m2d.txt",1,"number metals","prozent",1)
	
	#fig.add_subplot(2,2,4)
	#analyzer.plot_curve("mdis.txt",1,"number metals","prozent",1)
	plt.show()
