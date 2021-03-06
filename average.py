# average
import os
from analyzer import McAnalyzer
import platform
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from count_mol_bond import cal_bond_num, cluster
from PIL import Image

#filenumbers = [15,16,18,27,28,29,30,31,32,33,34,35,36,42,43,44,45]
filenumbers = [15,18,27,28,29,30,42,43,44,45]
temp_path = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/BDS286/results"

#totalenergy = np.zeros((9,13))
#totalenergy_av = np.zeros((9,13))
#mdense = np.zeros((9,13))
#m1d = np.zeros((9,13))
#m2d = np.zeros((9,13))
#mdis = np.zeros((9,13


# calculate the results
def totrun(temp_path,filenumbers):
    for temp_file in filenumbers:
        az = McAnalyzer('None')
        full_path = temp_path + str(temp_file)
        #print full_path
        az.set_path(full_path)
        print az.path
        az.load_logfile()
        az.run(1)


#totrun()
class McAverager:
    def __init__(self,wd_path,analyzer):
        self.path = wd_path
        self.analyzer = analyzer

    def set_path(self, new_path):
        self.path = new_path

    def set_filenumbers(self,filenumbers):
        self.filenumbers = filenumbers

    def init_range(self):
        fulldir = self.path + str(self.filenumbers[0])
        print fulldir
        os.chdir(fulldir)
        self.mdense = np.loadtxt("mdense.txt", delimiter=',')
        self.m1d = np.loadtxt("m1d.txt", delimiter=',')
        self.m2d = np.loadtxt("m2d.txt", delimiter=',')
        self.mdis = np.loadtxt("mdis.txt", delimiter=',')
        self.mtotal = np.loadtxt("mtotal.txt", delimiter=',')
        self.totalenergy = pd.Panel({str(self.filenumbers[0]):pd.DataFrame(np.loadtxt("totalenergy.txt",delimiter=','))})
        self.totalenergy_av = pd.Panel({str(self.filenumbers[0]):pd.DataFrame(np.loadtxt("totalenergy_av.txt",delimiter=','))})
        self.cbond_num = pd.Panel({str(self.filenumbers[0]):pd.DataFrame(np.loadtxt("cbond_num.txt",delimiter=','))})
        self.cbond_num_avt = pd.Panel({str(self.filenumbers[0]):pd.DataFrame(np.loadtxt("cbond_num_avt.txt",delimiter=','))})
        self.metal_ind = (self.mdense.shape)[0]
        self.energy_ind = (self.mdense.shape)[1]
        self.mdense_percent = self.mdense/self.mtotal
        self.m1d_percent = self.m1d/self.mtotal
        self.m2d_percent = self.m2d/self.mtotal

    def averg(self):
        for filenum in self.filenumbers[1:]:
            fulldir = self.path + str(filenum)
            print fulldir
            os.chdir(fulldir)
            temp_total = np.loadtxt("mtotal.txt", delimiter=',')
            temp_mdense = np.loadtxt("mdense.txt", delimiter=',')
            temp_m1d = np.loadtxt("m1d.txt", delimiter=',')
            temp_m2d = np.loadtxt("m2d.txt", delimiter=',')
            temp_mdis = np.loadtxt("mdis.txt", delimiter=',')
            temp_totalenergy = np.loadtxt("totalenergy.txt",delimiter=',')
            temp_totalenergy_av = np.loadtxt("totalenergy_av.txt",delimiter=',')
            temp_cbond_num = np.loadtxt("cbond_num.txt",delimiter=',')
            temp_cbond_num_avt = np.loadtxt("cbond_num_avt.txt",delimiter=',')

            self.mdense = self.mdense + temp_mdense
            self.m1d = self.m1d + temp_m1d
            self.m2d = self.m2d + temp_m2d
            self.mdis = self.mdis + temp_mdis
            self.totalenergy[str(filenum)] = pd.DataFrame(temp_totalenergy)
            self.totalenergy_av[str(filenum)] = pd.DataFrame(temp_totalenergy_av)
            self.totalenergy[str(filenum)] = pd.DataFrame(temp_totalenergy)
            self.totalenergy_av[str(filenum)] = pd.DataFrame(temp_totalenergy_av)
            self.cbond_num[str(filenum)] = pd.DataFrame(temp_cbond_num)
            self.cbond_num_avt[str(filenum)] = pd.DataFrame(temp_cbond_num_avt)
            #self.totalenergy = self.totalenergy + temp_totalenergy
            #self.totalenergy_av = self.totalenergy_av + temp_totalenergy_av
            self.mdense_percent = self.mdense_percent + temp_mdense/temp_total
            self.m1d_percent = self.m1d_percent + temp_m1d/temp_total
            self.m2d_percent = self.m2d_percent + temp_m2d/temp_total
            print str(filenum) + "is done"
        print "average is done!"
        #self.mdense = self.mdense/self.filenumbers
        os.chdir(self.path+str(self.filenumbers[0]))
        np.savetxt("mdense_final.txt",self.mdense,fmt='%0.2e',delimiter=',')
        np.savetxt("m1d_final.txt",self.m1d,fmt='%0.2e',delimiter=',')
        np.savetxt("m2d_final.txt",self.m2d,fmt='%0.2e',delimiter=',')
        np.savetxt("mdis_final.txt",self.mdis,fmt='%0.2e',delimiter=',')
        #np.savetxt("totalenergy_final.txt",self.totalenergy,fmt='0.2e',delimiter=',')
        #np.savetxt("totalenergy_av_final.txt",self.totalenergy_av,fmt='0.2e',delimiter=',')
        np.savetxt("mdense_p_final.txt",self.mdense_percent,fmt='%0.2e',delimiter=',')
        np.savetxt("m1d_p_final.txt",self.m1d_percent,fmt='%0.2e',delimiter=',')
        np.savetxt("m2d_p_final.txt",self.m2d_percent,fmt='%0.2e',delimiter=',')
        self.totalenergy.to_pickle("totalenergy_final.pkl")
        self.totalenergy_av.to_pickle("totalenergy_av_final.pkl")
        self.cbond_num.to_pickle("cbond_num_final.pkl")
        self.cbond_num_avt.to_pickle("cbond_num_avt_final.pkl")

    def process(self):
        try:
            # load the data
            os.chdir(self.path+str(self.filenumbers[0]))
            self.total = np.loadtxt("mtotal.txt", delimiter=',')
            self.mdense = np.loadtxt("mdense.txt", delimiter=',')
            self.m1d = np.loadtxt("m1d.txt", delimiter=',')
            self._m2d = np.loadtxt("m2d.txt", delimiter=',')
            self._mdis = np.loadtxt("mdis.txt", delimiter=',')
            self.metal_ind = (self.mdense.shape)[0]
            self.energy_ind = (self.mdense.shape)[1]
            #temp_totalenergy = np.loadtxt("totalenergy.txt",delimiter=',')
            #temp_totalenergy_av = np.loadtxt("totalenergy_av.txt",delimiter=',')
            self.totalenergy = pd.read_pickle("totalenergy_final.pkl")
            self.totalenergy_av = pd.read_pickle("totalenergy_av_final.pkl")
            self.cbond_num = pd.read_pickle("cbond_num_final.pkl")
            self.cbond_num_avt = pd.read_pickle("cbond_num_avt_final.pkl")
        except IOError as e:
            print "I/O error({0}):{1}".format(e.errno, e.strerror)
        else:
            print "Loading data done"
        # Begin to process the data
        # First <E^2>
        #E2_mean = np.zeros((self.metal_ind,self.energy_ind))
        #E_mean = np.zeros((self.metal_ind,self.energy_ind))
        #for i in self.filenumbers:
        #    temp_m1d = self.totalenergy[str(i)]
        # Energy part #########################################
        E2_mean = (self.totalenergy.pow(2)).sum(0)/self.totalenergy.shape[0]
        E_mean2 = ((self.totalenergy.sum(0))/self.totalenergy.shape[0]).pow(2)
        np.savetxt("E2_mean.txt",E2_mean,fmt='%0.2e',delimiter=',')
        np.savetxt("E_mean2.txt",E_mean2,fmt='%0.2e',delimiter=',')
        E_dif = E2_mean - E_mean2
        E_dif_av_n = np.zeros(E_dif.shape)
        #for i in range(E_dif.shape[0]):
        #    E_dif_av_n[i,:] = E_dif[i,:]/(300 + i*25)
        np.savetxt("E_dif.txt",E_dif,fmt='%0.2e',delimiter=',')
        np.savetxt("E_dif_av_n.txt",E_dif_av_n,fmt='%0.2e',delimiter=',')
        E2_mean_av = (self.totalenergy_av.pow(2)).sum(0)/self.totalenergy_av.shape[0]
        E_mean2_av = ((self.totalenergy_av.sum(0))/self.totalenergy_av.shape[0]).pow(2)
        np.savetxt("E2_mean_av.txt",E2_mean_av,fmt='%0.2e',delimiter=',')
        np.savetxt("E_mean2_av.txt",E_mean2_av,fmt='%0.2e',delimiter=',')
        E_dif_av_n2 = E2_mean_av - E_mean2_av
        np.savetxt("E_dif_av_n2.txt",E_dif_av_n2,fmt='%0.2e',delimiter=',')
        # coordination bonding part #############################
        C2_mean_avt = (self.cbond_num_avt.pow(2)).sum(0)/self.cbond_num_avt.shape[0]
        C_mean2_avt = ((self.cbond_num_avt.sum(0))/self.cbond_num_avt.shape[0]).pow(2)
        np.savetxt("C2_mean_avt.txt",C2_mean_avt,fmt='%0.2e',delimiter=',')
        np.savetxt("C_mean2_avt.txt",C_mean2_avt,fmt='%0.2e',delimiter=',')
        C_dif_avt = C2_mean_avt - C_mean2_avt
        np.savetxt("C_dif_avt.txt",C_dif_avt,fmt='%0.2e',delimiter=',')
        print "Process done"



if __name__ == '__main__':
    #totrun(temp_path,[42,43,44,45])
    az = McAnalyzer('None')
    averager = McAverager(temp_path,az)
    averager.set_filenumbers(filenumbers)
    averager.init_range()
    averager.averg()
    averager.process()
