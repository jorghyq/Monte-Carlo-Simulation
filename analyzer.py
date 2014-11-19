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

# go to the working directory
dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results13"
os.chdir(dname)
# load the data
files = os.listdir(dname)

def sprint(txt, inp):
	print txt + " : " + str(inp)


class McAnalyzer:
	def __init__(self, wd_path):
		self.path = wd_path
		
	def set_initial(self,nmet_init,nmet_step,cenergy_init,cernegy_step):
		self.nmet_init = nmet_init
		self.nmet_step = nmet_step
		self.cenergy_init = cenergy_init
		self.cernegy_step = cernegy_step
		
	
	def load_txt(self,txt_name):
		####### information from the name #######
		namedata = txt_name[0:-4].strip().split('-')
		self.latt_len = int(namedata[1])
		self.num_metal = int(namedata[3])
		self.nmetal_ind = (self.num_metal - self.nmet_init)/self.nmet_step
		self.cenergy = float(namedata[4])
		self.cenergy_ind = (self.cenergy - self.cenergy_init)/self.cernegy_step
		####### information from the header #######
		f = open(txt_name, 'r')
		headdata = f.readline().strip().split(',')
		f.close()
		self.cbond_num = float(headdata[0])
		self.cbond_num_av = self.cbond_num/float(self.num_metal)	
		self.vbond_num = float(headdata[1])
		self.total_energy = float(headdata[2])
		####### information from the matrix #######
		self.lattice = np.loadtxt(txt_name, delimiter=',',skiprows=1)
		self.lattice = self.lattice[0:self.latt_len,0:]
		
	
	def bond_num(self):
		temp = cal_bond_num(self.lattice) # 1*5 vector
		return temp
	
	def clustering(self,mode):
		mols,count = cluster(mode,self.lattice) # mols: list, count: total number
		return mols,count 

#class Summarizer:
#	def __init__(self):

nmetal_initial = 50
cenergy_initial = 1
nmetal_final = 550
cenergy_final = 19
nmetal_step = 25
cenergy_step = 2

ind_metal = (nmetal_final - nmetal_initial)/nmetal_step + 1
ind_cenergy = (cenergy_final - cenergy_initial)/cenergy_step + 1

analyzer = McAnalyzer(dname)
analyzer.set_initial(nmetal_initial,nmetal_step,cenergy_initial,cenergy_step)
#sprint("nmetal_intial",nmetal_initial)
#sprint("nmetal_final",nmetal_final)
#sprint("cenergy_initial",cenergy_initial)
#sprint("cenergy_final",cenergy_final)
#sprint("ind_metal",ind_metal)
#sprint("ind_cenergy",ind_cenergy)

#######################################
# dmol_num: ind_metal * ind_cenergy
#   1) dense
#   2) 1d
#   3) 2d
#	4) disordered
########################################	


#mdense = np.zeros((ind_metal,ind_cenergy))
#m1d = np.zeros((ind_metal,ind_cenergy))
#m2d = np.zeros((ind_metal,ind_cenergy))
#mdis = np.zeros((ind_metal,ind_cenergy))
totalenergy = np.zeros((ind_metal,ind_cenergy))
cbond_num = np.zeros((ind_metal,ind_cenergy))
vbond_num = np.zeros((ind_metal,ind_cenergy))
cbond_num_av = np.zeros((ind_metal,ind_cenergy))

for line in files:
	if line[0] == '1':
		analyzer.load_txt(line)
		totalenergy[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.total_energy
		cbond_num[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num
		vbond_num[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.vbond_num
		cbond_num_av[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num_av
		
np.savetxt("totalenergy.txt",totalenergy,delimiter=',')	
np.savetxt("cbond_num.txt",cbond_num,delimiter=',')
np.savetxt("vbond_num.txt",vbond_num,delimiter=',')
np.savetxt("cbond_num_av.txt",cbond_num_av,delimiter=',')
#print mdense.shape
#for line in files:
	#if line[0] == '1':
		#analyzer.load_txt(line)
		#mols0, count0 = analyzer.clustering(0)
		#mdense[analyzer.nmetal_ind][analyzer.cenergy_ind] = count0
		#mols1, count1 = analyzer.clustering(1)
		#m1d[analyzer.nmetal_ind][analyzer.cenergy_ind] = count1
		#mols2, count2 = analyzer.clustering(2)
		#m2d[analyzer.nmetal_ind][analyzer.cenergy_ind] = count2
		#newlist = []
		#count3 = 0
		#mol_list = range(200)
		#for i in range(len(mols0)):
			#for j in range(len(mols0[i])):
				#if mols0[i][j] not in newlist:
					#newlist.append(mols0[i][j])
		#for i in range(len(mols1)):
			#for j in range(len(mols1[i])):
				#if mols1[i][j] not in newlist:
					#newlist.append(mols1[i][j])
		#for i in range(len(mols2)):
			#for j in range(len(mols2[i])):
				#if mols2[i][j] not in newlist:
					#newlist.append(mols2[i][j])
		#for i in range(200):
			#if i not in newlist:
				#count3 = count3 + 1
		#mdis[analyzer.nmetal_ind][analyzer.cenergy_ind] = count3
		#sprint("count3",count3)
		#print line + " is done..."


#np.savetxt("mdense.txt",mdense, delimiter=',')	
#np.savetxt("m1d.txt",m1d, delimiter=',')
#np.savetxt("m2d.txt",m2d, delimiter=',')
#np.savetxt("mdis.txt",mdis, delimiter=',')
#mdense = np.loadtxt("mdense.txt", delimiter=',')
#m1d = np.loadtxt("m1d.txt", delimiter=',')
#m2d = np.loadtxt("m2d.txt", delimiter=',')
#mdis = np.loadtxt("mdis.txt", delimiter=',')
#mtotal = mdense + m1d + m2d + mdis

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
fig.add_subplot(2,2,1)
for i in range(0,totalenergy.shape[1]):
	plt.plot(xmetal,totalenergy[:,i],label = 'cenerg = %d' % (i*2+1))
	plt.text(xmetal[int(ind_metal/2)],totalenergy[int(ind_metal/2),i], str(i*2+1),fontsize=8)
plt.legend(loc=0,fontsize = 8)
plt.xlabel('number of metals')
plt.ylabel('total energy')


fig.add_subplot(2,2,2)
for i in range(0,cbond_num.shape[1]):
	plt.plot(xmetal,cbond_num[:,i],label = 'cenerg = %d' % (i*2+1))
	plt.text(xmetal[int(ind_metal/2)],cbond_num[int(ind_metal/2),i], str(i*2+1),fontsize=8)
plt.legend(loc=0,fontsize = 8)
plt.xlabel('number of metals')
plt.ylabel('number of coordination bond')

fig.add_subplot(2,2,3)
for i in range(0,vbond_num.shape[1]):
	plt.plot(xmetal,vbond_num[:,i],label = 'cenerg = %d' % (i*2+1))
	plt.text(xmetal[int(ind_metal/2)],vbond_num[int(ind_metal/2),i], str(i*2+1),fontsize=8)
plt.legend(loc=0,fontsize = 8)
plt.xlabel('number of metals')
plt.ylabel('number of vdW bond')

fig.add_subplot(2,2,4)
for i in range(0,cbond_num_av.shape[1]):
	plt.plot(xmetal,cbond_num_av[:,i],label = 'cenerg = %d' % (i*2+1))
	plt.text(xmetal[int(ind_metal/2)],cbond_num_av[int(ind_metal/2),i], str(i*2+1),fontsize=8)
plt.legend(loc=0,fontsize = 8)
plt.xlabel('number of metals')
plt.ylabel('ratio between coordiation bond/number of metals')
plt.show()
