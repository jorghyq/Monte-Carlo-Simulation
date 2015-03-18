import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster
from PIL import Image

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results37"
os.chdir(dname)

cenergy_20 =  np.loadtxt("cenergy_20.txt", delimiter=',')
cenergy_10 =  np.loadtxt("cenergy_10.txt", delimiter=',')
cenergy_6 =  np.loadtxt("cenergy_6.txt", delimiter=',')
cenergy_5 =  np.loadtxt("cenergy_5.txt", delimiter=',')
cenergy_4 =  np.loadtxt("cenergy_4.txt", delimiter=',')
cenergy_3 =  np.loadtxt("cenergy_3.txt", delimiter=',')
#metal_650 =  np.loadtxt("metal_650.txt", delimiter=',')

mdense = np.zeros((6,9))
m1d = np.zeros((6,9))
m2d = np.zeros((6,9))
mdis = np.zeros((6,9))
i = 0
filename = "cenergy_20.txt"
print filename 
tmp = np.loadtxt(filename, delimiter=',')
mdense[i,:] = tmp[1,:]
m1d[i,:] = tmp[2,:]
m2d[i,:] = tmp[3,:]
mdis[i,:] = tmp[4,:]
i = i + 1
filename = "cenergy_10.txt"
print filename 
tmp = np.loadtxt(filename, delimiter=',')
mdense[i,:] = tmp[1,:]
m1d[i,:] = tmp[2,:]
m2d[i,:] = tmp[3,:]
mdis[i,:] = tmp[4,:]
i = i + 1
filename = "cenergy_6.txt"
print filename 
tmp = np.loadtxt(filename, delimiter=',')
mdense[i,:] = tmp[1,:]
m1d[i,:] = tmp[2,:]
m2d[i,:] = tmp[3,:]
mdis[i,:] = tmp[4,:]
i = i + 1
filename = "cenergy_5.txt"
print filename 
tmp = np.loadtxt(filename, delimiter=',')
mdense[i,:] = tmp[1,:]
m1d[i,:] = tmp[2,:]
m2d[i,:] = tmp[3,:]
mdis[i,:] = tmp[4,:]
i = i + 1
filename = "cenergy_4.txt"
print filename 
tmp = np.loadtxt(filename, delimiter=',')
mdense[i,:] = tmp[1,:]
m1d[i,:] = tmp[2,:]
m2d[i,:] = tmp[3,:]
mdis[i,:] = tmp[4,:]
i = i + 1
filename = "cenergy_3.txt"
print filename 
tmp = np.loadtxt(filename, delimiter=',')
mdense[i,:] = tmp[1,:]
m1d[i,:] = tmp[2,:]
m2d[i,:] = tmp[3,:]
mdis[i,:] = tmp[4,:]
	
np.savetxt("mdense.txt",mdense, delimiter=',')	
np.savetxt("m1d.txt",m1d, delimiter=',')
np.savetxt("m2d.txt",m2d, delimiter=',')
np.savetxt("mdis.txt",mdis, delimiter=',')

mtotal = mdense + m1d + m2d + mdis
np.savetxt("mtotal.txt",mtotal, delimiter=',')
print "Done, data are saved!"
	


def phase_diagram(updown,leftright,xlab,ylab):
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
	plt.imshow(im,extent=[0.0,1,0.05,0.3],aspect="auto")
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.show()

phase_diagram(1,1,"Ev/Ec","kT/Ec")
