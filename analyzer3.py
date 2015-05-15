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
		
	def load_txt(self,txt_name):
		self.current_mctext = McText(txt_name)
	
	def run(self,num_metal):
		lines = os.listdir(self.path)
		os.chdir(dname)
		self.result = np.zeros((5,7))
		self.result2 = np.zeros((5,7))
		p = 0
		for line in lines:
			if line[0] == '1' and line[-4:] == '.txt':
				
				self.load_txt(line)
				if self.current_mctext.num_metal == num_metal:
					print line
					print self.current_mctext.num_metal
					self.result[0][p] = self.current_mctext.venergy
					mols0, count0 = self.current_mctext.clustering(0)
					self.result[1][p] = count0
					mols1, count1 = self.current_mctext.clustering(1)
					self.result[2][p] = count1
					mols2, count2 = self.current_mctext.clustering(3)
					self.result[3][p] = count2
					newlist = []
					count3 = 0
					mol_list = range(self.current_mctext.num_mol)
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
					for i in range(self.current_mctext.num_mol):
						if i not in newlist:
							count3 = count3 + 1
					self.result[4][p] = count3
					p = p + 1
		index = np.argsort(-self.result[0,:])
		for i in range(0,5):
			self.result2[i,:] = self.result[i,index]
		np.savetxt("metal_%d.txt" % num_metal,self.result2,delimiter=',')
				
			
class McText:
	def __init__(self,txt_name):
		self.txt_name = txt_name
		temp, filename = os.path.split(txt_name)
		namedata = filename[0:-4].strip().split('_')
		self.latt_len = int(namedata[1])
		self.num_metal = int(namedata[3])
		self.num_mol = int(namedata[2])
		#self.nmetal_ind = (self.num_metal - self.nmet_init)/self.nmet_step
		self.cenergy = float(namedata[4])
		self.venergy = float(namedata[5])
		#self.cenergy_ind = (self.cenergy - self.cenergy_init)/self.cenergy_step
		#self.venergy_ind = (self.venergy - self.venergy_init)/self.venergy_step
		####### information from the header #######
		f = open(txt_name, 'r')
		headdata = f.readline().strip().split(',')
		f.close()
		#self.cbond_num = float(headdata[0])
		#self.cbond_num_avt = self.cbond_num/float(self.num_metal + self.num_mol)	
		#self.vbond_num = float(headdata[1])
		#self.total_energy = float(headdata[2])
		#self.energy_av = self.total_energy/float(self.num_metal + self.num_mol)
		####### information from the matrix #######
		self.lattice = np.loadtxt(txt_name, delimiter=',',skiprows=1)
		self.lattice = self.lattice[0:self.latt_len,0:]
		#sprint("num_metal",self.num_metal)

	def bond_num(self):
		temp = cal_bond_num(self.lattice) # 1*5 vector
		return temp
	
	def clustering(self,mode):
		mols,count = cluster(mode,self.lattice) # mols: list, count: total number
		return mols,count
				

if __name__ == "__main__":
	# go to the working directory
	dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results41"
	analyzer = McAnalyzer(dname)
	for i in range(0,360,50):
		analyzer.run(i)
	#analyzer.run(0)
	#print analyzer.result
	#print analyzer.result2
	#fig = plt.figure()
	#plt.show()
