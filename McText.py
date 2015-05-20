# Python Monte Carlo Simulation Classes
# This file incluees: McText, McLog
#
import os
import numpy as np
# import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster


class McText:
    def __init__(self):
        self.mol_type = 0  # default BDS286 = 0, BDS285 = 1
        pass

    def loadtxt(self, full_path):
        # TODO: Automatic identify the species
        self.full_path = full_path
        self.path, self.txt_name = os.path.split(self.full_path)
        # information from the file name ######
        namedata = self.txt_name[0:-4].strip().split('_')
        self.latt_len = int(namedata[1])
        self.num_metal = int(namedata[3])
        self.num_mol = int(namedata[2])
        self.cenergy = float(namedata[4])
        self.venergy = float(namedata[5])
        # information from the header #######
        f = open(self.txt_name, 'r')
        headdata = f.readline().strip().split(',')
        f.close()
        self.cbond_num = float(headdata[0])
        self.cbond_num_avt = self.cbond_num/float(self.num_metal + self.num_mol)
        self.vbond_num = float(headdata[1])
        self.total_energy = float(headdata[2])
        self.energy_av = self.total_energy/float(self.num_metal + self.num_mol)
        # information from the matrix #######
        self.lattice = np.loadtxt(self.txt_name, delimiter=',', skiprows=1)
        self.lattice = self.lattice[0:self.latt_len, 0:]

    def bond_num(self):
        temp = cal_bond_num(self.lattice)
        return temp

    def mol_cluster(self):
        if self.mol_type == 0:
            mols0, count0 = cluster(0, self.lattice)
            self.mdense = count0
            mols1, count1 = cluster(1, self.lattice)
            self.m1d = count1
            mols2, count2 = cluster(2, self.lattice)
            self.m2d = count2
            newlist = []
            count3 = 0
            mol_list = range(self.num_mol)
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
            for i in range(self.num_mol):
                if i not in newlist:
                    count3 = count3 + 1
            self.mdis = count3
        elif self.mol_type == 1:
            mols0, count0 = cluster(0, self.lattice)
            self.mdense= count0
            mols1, count1 = cluster(1, self.lattice)
            self.m1d = count1
            mols2, count2 = cluster(3, self.lattice)
            self.m2d = count2
            newlist = []
            count3 = 0
            mol_list = range(self.num_mol)
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
            for i in range(self.num_mol):
                if i not in newlist:
                    count3 = count3 + 1
            self.mdis = count3
        mtotal = mdense + m1d + m2d + mdis
        self.mdense_p = self.mdense / self.mtotal
        self.m1d_p = self.m1d / self.mtotal
        self.m2d_p = self.m2d / self.mtotal
        self.mdis_p = self.mdis / self.mtotal

    def energy_cal(self):
        pass

    def process(self):
        # Calculate the molecules belonging to each catogory
        pass

class McLog:
    def __init__(self,path):
        self.path = path
        self.mol_type = 0

    def init_dir(self,path): # Init the logfile from the files inside current directory
        #self.total_run
        pass

    def set_param(self,param): # param=[total_run,latt_len,num_mol,nmet_num,cenergy,venergy]
        self.total_run = param[0]
        self.latt_len = param[1]
        self.num_mol = param[2]
        self.nmet_init = param[3]
        self.nmet_step = param[4]
        self.nmet_max = param[5]
        self.cenergy_init = param[6]
        self.cenergy_step = param[7]
        self.cenergy_max = param[8]
        self.venergy_init = param[9]
        self.venergy_step = param[10]
        self.venergy_max = param[11]

    def load_log(self):
        fname = os.path.join(self.path,'logfile.txt')
        #print fname
        if os.path.isfile(fname):
            with open(fname, 'r') as f:
                print "load logfile!"
                #counter = 0
                for line in f:
                    tline = line.split(':')
                    #print tline
                    if len(tline) > 1:
                        if tline[0].strip() == 'total_run':
                            self.total_run = float(tline[1].strip())
                        elif tline[0].strip() == 'mol_type':
                            self.mol_type = int(tline[1].strip())
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
                            print "logfile loaded!"
                            self.ind_metal = (self.nmet_max - self.nmet_init)/self.nmet_step + 1
                            self.ind_cenergy = (self.cenergy_max - self.cenergy_init)/self.cenergy_step + 1
                            self.ind_venergy = (self.venergy_max - self.venergy_init)/self.venergy_step + 1
                        else:
                            print "no logfile!!!!!"

    def save_log(self, comment):
        with open(os.path.join(self.path,'logfile.txt'), 'w') as f:
            f.write('This is the logfile for the following settings.\n')
            f.write(comment + '.\n')
            f.write('mol_type: '+ str(self.mol_type))
            f.write('total_run: ' + str(self.total_run) + '\n')
            f.write('latt_len: ' + str(self.latt_len) + '\n')
            f.write('num_mol: ' + str(self.num_mol) + '\n')
            f.write('nmet_init: ' + str(self.nmet_init) + '\n')
            f.write('nmet_step: ' + str(self.nmet_step) + '\n')
            f.write('nmet_max: ' + str(self.nmet_max) + '\n')
            f.write('cenergy_init: ' + str(self.cenergy_init) + '\n')
            f.write('cenergy_step: ' + str(self.cenergy_step) + '\n')
            f.write('cenergy_max: ' + str(self.cenergy_max) + '\n')
            f.write('venergy_init: ' + str(self.venergy_init) + '\n')
            f.write('venergy_step: ' + str(self.venergy_step) + '\n')
            f.write('venergy_max: ' + str(self.venergy_max) + '\n')
            print "logfile saved to" + os.path.join(self.path,'logfile.txt')""
