# Python file to analyze the lattice obtained from the simulation with the form specified below
# the center of the molecules : 3
# legs : 2
# metals : 1

import os
import platform
import sys
import numpy as np
import matplotlib.pyplot as plt
#from count_mol_bond import cal_bond_num, cluster
from McText import McText, McLog
from PIL import Image


def sprint(txt, inp):
    print txt + " : " + str(inp)


class McAnalyzer:
    def __init__(self, wd_path):
        self.path = wd_path
        self.mct = McText()
        self.mcl = McLog()
        self.mol_type = 0
        self.mct.mol_type = self.mol_type

    def set_path(self, new_path):
        self.path = new_path

    def set_initial(self, nmet_init, nmet_step, nmet_max, \
                    cenergy_init, cenergy_step, cenergy_max, \
                    venergy_init, venergy_step, venergy_max):
        self.nmet_init = nmet_init
        self.nmet_step = nmet_step
        self.nmet_max = nmet_max
        self.cenergy_init = cenergy_init
        self.cenergy_step = cenergy_step
        self.cenergy_max = cenergy_max
        self.venergy_init = venergy_init
        self.venergy_step = venergy_step
        self.venergy_max = venergy_max
        self.ind_metal = (self.nmet_max - self.nmet_init)/self.nmet_step + 1
        self.ind_cenergy = (self.cenergy_max - self.cenergy_init)/self.cenergy_step + 1
        self.ind_venergy = (self.venergy_max - self.venergy_init)/self.venergy_step + 1
        if self.cenergy_init == self.cenergy_max:
            index =  range(self.nmetal_init, self.nmetal_max, self.nmetal_max+self.nmetal_step)
            columns = range(self.venergy_init, self.venergy_max, self.venergy_max+self.venergy_step)
            temp = pd.DataFrame(index=index, columns=columns)
            temp = temp.fillna(0)
        # elif self.venergy_init == self.venergy_max:
        #    index =
        #    column =
        self.data = pd.Panel({'proto':temp})

    def run(self,mode):
        files = os.listdir(self.path)
        os.chdir(self.path)
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
                            #if os.path.exists(os.path.join(self.path,"mdense.txt")):
                            #	mdis = np.loadtxt("mdis.txt", delimiter=',')
                            print "DATA ALREADY EXIST!"
                            return
            print "Begin to process the data..."
            for line in files:
                if line[0] == '1' and line[-4:] == '.txt':
                    self.mct.load_txt(line)
                    # determine the index of the current file
                    cenergy_ind = (self.mct.cenergy - self.cenergy_init)/self.cenergy_step
                    venergy_ind = (self.mct.venergy - self.venergy_init)/self.venergy_step
                    nmetal_ind = (self.mct.num_metal - self.nmet_init)/self.nmet_step
                    totalenergy[self.nmetal_ind][self.venergy_ind] = self.total_energy
                    cbond_num[self.nmetal_ind][self.venergy_ind] = self.cbond_num
                    vbond_num[self.nmetal_ind][self.venergy_ind] = self.vbond_num
                    #cbond_num_av[analyzer.nmetal_ind][analyzer.cenergy_ind] = analyzer.cbond_num_av
                    totalenergy_av[self.nmetal_ind][self.venergy_ind] = self.energy_av
                    cbond_num_avt[self.nmetal_ind][self.venergy_ind] = self.cbond_num_avt
                    if mode == 1:
                        self.mct.
            np.savetxt("totalenergy.txt",totalenergy,delimiter=',',fmt='%0.4f')
            np.savetxt("cbond_num.txt",cbond_num,delimiter=',',fmt='%0.4f')
            np.savetxt("vbond_num.txt",vbond_num,delimiter=',',fmt='%0.4f')
            np.savetxt("cbond_num_av.txt",cbond_num_av,delimiter=',',fmt='%0.4f')
            np.savetxt("totalenergy_av.txt",totalenergy_av,delimiter=',',fmt='%0.4f')
            np.savetxt("cbond_num_avt.txt",cbond_num_avt,delimiter=',',fmt='%0.4f')
            np.savetxt("mdense.txt",mdense, delimiter=',',fmt='%0.4f')
            np.savetxt("m1d.txt",m1d, delimiter=',',fmt='%0.4f')
            np.savetxt("m2d.txt",m2d, delimiter=',',fmt='%0.4f')
            np.savetxt("mdis.txt",mdis, delimiter=',',fmt='%0.4f')
            np.savetxt("mtotal.txt",mtotal, delimiter=',',fmt='%0.4f')
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
                plt.imshow(im,extent=[0.125,1.125,0/self.num_mol,600/self.num_mol],aspect="auto")
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
    filenumber = [46]#,27,28,29,30,31,33]
    for filenum in filenumber:
        dname = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results"+str(filenum)
        print "file " + str(filenum) + " in processing"
        analyzer = McAnalyzer(dname)
        analyzer.load_logfile()
        analyzer.run(2)
        #fig = plt.figure()
        #analyzer.phase_diagram(1,0,"Ev/Ec","metal/molecule")
        #fig.add_subplot(2,2,1)
        #analyzer.plot_curve("mdense.txt",1,"number metals","prozent",1)

    #fig.add_subplot(2,2,2)
    #analyzer.plot_curve("m1d.txt",1,"number metals","prozent",1)

    #fig.add_subplot(2,2,3)
    #analyzer.plot_curve("m2d.txt",1,"number metals","prozent",1)

    #fig.add_subplot(2,2,4)
    #analyzer.plot_curve("mdis.txt",1,"number metals","prozent",1)
    plt.show()
