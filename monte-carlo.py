# Monte Carlo
import math
import numpy as np
import scipy as sp
import random as rd
import matplotlib.pyplot as plt

# define some unchanged variables
global LATTICE_CONSTANT
LATTICE_CONSTANT = 0.2
print LATTICE_CONSTANT



########################################################################
# define the lattice class
class Lattice:
	global LATTICE_CONSTANT
	def __init__(self, height, width, angle):
		self.height = height
		self.width = width
		# define the angle of the unit cell, assume all the side of the unit cell is the same
		self.angle = angle
		# define the lattice points 
		self.lattx = np.zeros((height, width))
		self.latty = np.zeros((height, width))
		tempx = 0
		tempy = 0
		for i in range(0, self.height):
			tempx = LATTICE_CONSTANT * math.cos(self.angle) * i
			for j in range(0, self.width):
				self.lattx[i, j] = tempx
				self.latty[i, j] = tempy
				tempx = tempx + LATTICE_CONSTANT
			tempy = tempy - LATTICE_CONSTANT * math.sin(self.angle)			
			#print tempx
		# define the state of each points, 0 for not ocupied, 1 for ocupied
		self.lattstate = np.zeros((height, width))
		
	def tell(self):
		print self.lattx
		print self.latty

	def display(self):
		plt.figure
		for i in range(0, self.height):
			plt.scatter(self.lattx[i,:], self.latty[i,:])
		xmax = np.amax(self.lattx)
		print xmax
		xmin = np.amin(self.lattx)
		ymax = np.amax(self.latty)
		ymin = np.amin(self.latty)
		print xmin
		print ymax
		print ymin
		#elem = [xmax, xmin, ymax, ymin]
		pmin = min(xmin, ymin)
		pmax = max(xmax, ymax)
		print pmin
		print pmax
		plt.xlim(pmin, pmax)
		plt.ylim(pmin, pmax)
		#plt.show()

	#def change(self,)
		
########################################################################
# define the molecular class and the metal class	


class Things():
	def __init__(self, coordinate):
		self.coordinate = coordinate

	def newplace(self, coordinate):
		self.coordinate = coordinate

class Metal(Things):
	def __init__(self, coordinate):
		Things.__init__(self, coordinate)

	def newplace(self, coordinate):
		Things.newplace(self, coordinate)
		
class Molecule(Things):
	global LATTICE_CONSTANT
	def __init__(self, coordinate, angle):
		Things.__init__(self, coordinate)
		# define the angle of the molecule
		self.angle = angle
		# define the nodes of the molecule
		#self.numofpoints = numberofpoints
		self.points = np.zeros((3,2))
		self.prefered = np.zeros((2,2))
		#print coordinate[0]
		tempx = coordinate[0] + LATTICE_CONSTANT * math.sin(angle)
		tempy = coordinate[1] + LATTICE_CONSTANT * math.cos(angle)
		print tempx
		print tempy
		self.points[0,:] = np.array([(coordinate[0] + LATTICE_CONSTANT * math.sin(angle)), (coordinate[1] + LATTICE_CONSTANT * math.cos(angle))])
		self.points[1,:] = self.coordinate
		self.points[2,:] = np.array([coordinate[0] - LATTICE_CONSTANT * math.sin(angle), coordinate[1] - LATTICE_CONSTANT * math.cos(angle)])
		#self.prefered[1,:] = 

	def newplace(self, coordinate, angle):
		Things.newplace(self, coordinate)
		self.points[0,:] = np.array([coordinate[0] + LATTICE_CONSTANT * math.sin(angle), coordinate[1] + LATTICE_CONSTANT * math.cos(angle)])
		self.points[1,:] = coordinate
		self.points[2,:] = np.array([coordinate[0] - LATTICE_CONSTANT * math.sin(angle), coordinate[1] - LATTICE_CONSTANT * math.cos(angle)])

	def display(self):
		plt.plot([self.points[0,0],self.points[2,0]],[self.points[0,1],self.points[2,1]],linewidth=5.0)
		#plt.show()

########################################################################
	
# define the experiment class

class Experiment:
	#lattice
	def __init__(self, lattice, angle, numMolecule, numMetal):
		self.lattice = lattice
		self.nMo = numMolecule
		self.nMe = numMetal
		self.totalRun = 1000
		self.molecules = np.ndarray((self.nMo-1,),dtype = np.object)
		self.metals = np.ndarray((self.nMe-1,),dtype = np.object)
		self.angleGroup = np.arange(0,360,angle)
	
	#def calculate(self, indexOfMolecule, indexOfMetal):
		
		
	def initialize(self):
		# put all the molecules and metals on the lattice
		for i in range(0, self.nMo-1):
			indX = rd.randint(0, self.lattice.width-1)
			indY = rd.randint(0, self.lattice.height-1)
			indAngle = rd.randint(0,self.angleGroup.size-1)
			tempCoordinate = np.array([self.lattice.lattx[indX,indY],self.lattice.latty[indX,indY]])
			self.molecules[i] = Molecule(np.array([self.lattice.lattx[indX,indY],self.lattice.latty[indX,indY]]),self.angleGroup[indAngle])
			#print tempCoordinate
			
			
		
########################################################################
# main program

# create a lattice		
lat = Lattice(30,30, 3.14/3)
exp = Experiment(lat, 60, 10, 10)
exp.initialize()
exp.lattice.display()
for i in range(0, exp.molecules.size):
	exp.molecules[i].display()
#exp.molecules[0] = Molecule([[1],[2]],60)
#print exp.molecules[0].points
# create an experiment
plt.show()

