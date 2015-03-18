import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster
from PIL import Image

dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results41"
os.chdir(dname)

metal_0 =  np.loadtxt("metal_0.txt", delimiter=',')
metal_50 =  np.loadtxt("metal_50.txt", delimiter=',')
metal_100 =  np.loadtxt("metal_100.txt", delimiter=',')
metal_150 =  np.loadtxt("metal_150.txt", delimiter=',')
metal_200 =  np.loadtxt("metal_200.txt", delimiter=',')
metal_250 =  np.loadtxt("metal_250.txt", delimiter=',')
metal_300 =  np.loadtxt("metal_300.txt", delimiter=',')
metal_350 =  np.loadtxt("metal_350.txt", delimiter=',')
#metal_400 =  np.loadtxt("metal_400.txt", delimiter=',')
#metal_450 =  np.loadtxt("metal_450.txt", delimiter=',')
#metal_500 =  np.loadtxt("metal_500.txt", delimiter=',')
#metal_550 =  np.loadtxt("metal_550.txt", delimiter=',')
#metal_600 =  np.loadtxt("metal_600.txt", delimiter=',')
#metal_650 =  np.loadtxt("metal_650.txt", delimiter=',')

mdense = np.zeros((8,7))
m1d = np.zeros((8,7))
m2d = np.zeros((8,7))
mdis = np.zeros((8,7))

for i in range(8):
	filename = "metal_%d.txt" % (i*50)
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
	plt.imshow(im,extent=[0.025,0.325,0/300,350/300],aspect="auto")
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.show()

phase_diagram(1,0,"kT/Ec","metal/molecule")
