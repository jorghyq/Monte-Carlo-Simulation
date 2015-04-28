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

filenumbers = [15,16,18,27,28,29,30,31,32,33,34,35,36,42,43,44,45]

temp_path = "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results"

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
        #self.analyzer.set_path(self.path)
        #self.analyzer.load_logifle()

        # read all the files
        #files = os.listdir(self.path)
        # select the files
        #for file in files:
        #    if file == "logfile.txt":



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
        self.metal_ind = (self.mdense.shape)[0]
        self.energy_ind = (self.mdense.shape)[1]
        self.mdense_percent = self.mdense/self.mtotal
        self.m1d_percent = self.m1d/self.mtotal
        self.m2d_percent = self.m2d/self.mtotal

    def averg(self):
        for filenum in self.filenumbers:
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
            self.mdense = self.mdense + temp_mdense
            self.m1d = self.m1d + temp_m1d
            self.m2d = self.m2d + temp_m2d
            self.mdis = self.mdis + temp_mdis
            self.totalenergy[str(filenum)] = pd.DataFrame(temp_totalenergy)
            self.totalenergy_av[str(filenum)] = pd.DataFrame(temp_totalenergy_av)
            #self.totalenergy = self.totalenergy + temp_totalenergy
            #self.totalenergy_av = self.totalenergy_av + temp_totalenergy_av
            self.mdense_percent = self.mdense_percent + temp_mdense/temp_total
            self.m1d_percent = self.m1d_percent + temp_m1d/temp_total
            self.m2d_percent = self.m2d_percent + temp_m2d/temp_total
        print "average is done!"
        self.mdense = self.mdense/self.filenumbers
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
        E2_mean = (self.totalenergy.pow(2)).sum(0)/self.totalenergy.shape[0]
        E_mean2 = ((self.totalenergy.sum(0))/self.totalenergy.shape[0]).pow(2)
        np.savetxt("E2_mean.txt",E2_mean,fmt='%0.2e',delimiter=',')
        np.savetxt("E_mean2.txt",E_mean2,fmt='%0.2e',delimiter=',')
        E_dif = E2_mean - E_mean2
        np.savetxt("E_dif.txt",E_dif,fmt='%0.2e',delimiter=',')
        print "Process done"



if __name__ == '__main__':
    totrun(temp_path,[16])
    #az = McAnalyzer('None')
    #averager = McAverager(temp_path,az)
    #averager.set_filenumbers(filenumbers)
    #averager.init_range()
    #averager.averg()
    #averager.process()
    #totrun(temp_path,[33])







