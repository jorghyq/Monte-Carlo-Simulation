# Simply convert all the text files to the image file
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
if len(sys.argv) == 2:
    ffn = int(sys.argv[1])

#dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results47"
dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results" + str(ffn)+'/'
print dname
os.chdir(dname)
files = os.listdir(dname)

for line in files:
    if (line[-4:] == '.txt' and line[0] == '1'):
        newline = dname+line[:-4]+'.png'
        print line
        print newline
        namedata = line[0:-4].strip().split('_')
    	latt_len = int(namedata[1])
    	print latt_len
        lattice = np.loadtxt(dname+line, skiprows = 1,delimiter=',')
        lattice = lattice[0:latt_len,:]
        print lattice.shape
        #plt.imsave(newline,lattice,[0,2])
        temp1, temp2 = np.where(lattice == 9) # mol1
        temp3, temp4 = np.where(lattice == 10) # mol2
        temp5, temp6 = np.where(lattice == 1)  # metal
        temp7, temp8 = np.where(lattice == 7)
        # To have a customer designed shape, one has to draw it by himself
        x = [-1.5,-0.5,-0.5,0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5]
        y = [0.5,0.5,1.5,1.5,0.5,0.5,-0.5,-0.5,-1.5,-1.5,-0.5,-0.5,0.5]
        x2 = [-3,-1,-1,1,1,3,3,1,1,-1,-1,-3,-3]
        y2 = [1,1,3,3,1,1,-1,-1,-3,-3,-1,-1,1]
        x3 = [-0.5,0.5,0.5,-0.5,-0.5]
        y3 = [0.5,0.5,-0.5,-0.5,-0.5]
        #x3 = x/2
        #y3 = y/2
        xy1 = list(zip(x,y))
        xy2 = list(zip(x2,y2))
        xy3 = list(zip(x3,y3))
        fig = plt.figure()
        ax = plt.axes()
        fig.add_axes(ax)

        ax.scatter(temp1,temp2,s = 120,c = "#0ACEF5",linewidth='0',marker = xy1)
        ax.scatter(temp3,temp4,s = 120,c = "b",linewidth='0',marker = xy1)
        ax.scatter(temp7,temp8,s = 20,c = "r",linewidth='0',marker = xy3)
        ax.scatter(temp5,temp6,s = 20,c = "#F78C00",linewidth='0',marker = "o")

        plt.xlim([-0.5, latt_len-0.5])
        plt.ylim([-0.5, latt_len-0.5])

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        #ax.text(-10,0,headdate[0])
        ax.set_aspect("equal")
        #plt.show()
        plt.savefig(newline,dpi=500,bbox_inches='tight',pad_inches=0)

