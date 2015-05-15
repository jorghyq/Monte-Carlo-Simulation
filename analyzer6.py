import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster
from PIL import Image

dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results54"
os.chdir(dname)
cenergy = [40,20,13,10,8.0,6.7,5.7,5,4.4,4,3.6,3.3]
cenergy_length = len(cenergy)
mdense = np.zeros((cenergy_length,9))
m1d = np.zeros((cenergy_length,9))
m2d = np.zeros((cenergy_length,9))
mdis = np.zeros((cenergy_length,9))

for i in range(cenergy_length):
	t_filename = 'cenergy_%0.1f.txt' % cenergy[i]
	tmp = np.loadtxt(t_filename, delimiter=',')
	mdense[i,:] = tmp[1,:]
	m1d[i,:] = tmp[2,:]
	m2d[i,:] = tmp[3,:]
	mdis[i,:] = tmp[4,:]
	print t_filename + " readed!"

np.savetxt("mdense.txt",mdense, delimiter=',',fmt='%0.4f')	
np.savetxt("m1d.txt",m1d, delimiter=',',fmt='%0.4f')
np.savetxt("m2d.txt",m2d, delimiter=',',fmt='%0.4f')
np.savetxt("mdis.txt",mdis, delimiter=',',fmt='%0.4f')

mtotal = mdense + m1d + m2d + mdis
np.savetxt("mtotal.txt",mtotal, delimiter=',',fmt='%0.4f')
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
	plt.imshow(im,extent=[0.0,1,0.025,0.3],aspect="auto")
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	out_name = 'BDS285_kT_vs_Eo.png'
	plt.savefig(out_name,dpi=500)
	plt.show()

phase_diagram(1,1,"Eo/Ec","kT/Ec")

