# Simply convert all the text files to the image file
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

s_cross = 120
s_rec  = 20

if len(sys.argv) == 2:
    ffn = int(sys.argv[1])
elif len(sys.argv) == 4:
    ffn = int(sys.argv[1])
    s_cross = int(sys.argv[2])
    s_rec = int(sys.argv[3])
#dname = "D:\Dropbox\Project\python\Monte-Carlo-Simulation\\results47"
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

with open('colors.json','r') as f:
    colors = json.load(f)

dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results" + str(ffn)+'/'
print dname
os.chdir(dname)
files = os.listdir(dname)

for line in files:
    if (line[-4:] == '.txt' and line[0] == '1' and line[-5] != 't'):
        newline = dname+line[:-4]+'.png'
        print line
        print newline
        namedata = line[0:-4].strip().split('_')
    	latt_len = int(namedata[1])
    	print latt_len
        lattice = np.loadtxt(dname+line, skiprows = 1,delimiter=',')
        lattice = lattice[0:latt_len,:]
        print lattice.shape
        fig = plt.figure()
        ax = plt.axes()
        fig.add_axes(ax)
        coor_x, coor_y = np.where(lattice == 1) # metal
        ax.scatter(coor_x,coor_y,s=s_rec,c=colors['1'],linewidth='0',marker = "o")

        for i in [11,12]:
            coor_x,coor_y = np.where(lattice == i)
            ax.scatter(coor_x,coor_y,s=s_cross,c=colors[str(i)],linewidth='0',marker=xy1)
        # plot endgroups
        for i in range(2,11):
            coor_x,coor_y = np.where(lattice == i)
            if str(i) in colors:
                ax.scatter(coor_x,coor_y,s=s_rec,c=colors[str(i)],linewidth='0',marker=xy3)
        plt.xlim([-0.5, latt_len-0.5])
        plt.ylim([-0.5, latt_len-0.5])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        #ax.text(-10,0,headdate[0])
        ax.set_aspect("equal")
        #plt.show()#plt.imsave(newline,lattice,[0,2])
        plt.savefig(newline,dpi=500,bbox_inches='tight',pad_inches=0)
        #plt.show()

