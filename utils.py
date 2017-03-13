"""
functionsto analysis simulation results
molecular center: 11, 12

"""
import numpy as np

m_center = [11, 12]
directions = np.array([[-1,0],[0,1],[1,0],[0,-1]])

def kwn(n, var):
    """ keep var within n """
    if var > n-1:
        temp = var - n
    elif var < 0:
        temp = var + n
    else:
        temp = var
    return temp

def get_neighbour(pos, lat_len, order):
    #print pos
    """ return the pos around the input pos """
    neighbours = np.tile(pos,(4,1)) + order * directions
    # limit the coordinates within [0, lat_len-1]
    x,y = np.where(neighbours > lat_len - 1)
    coor = zip(x,y)
    for item in coor:
        neighbours[tuple(item)] = neighbours[tuple(item)] - lat_len
    x,y = np.where(neighbours < 0)
    coor = zip(x,y)
    for item in coor:
        neighbours[tuple(item)] = neighbours[tuple(item)] + lat_len
    return neighbours

def num_mbond(lattice, end_group):
    """ return the number of metal-coordination bond formed by end_group """
    count = 0
    # get the coordinates of the molecules
    x1, y1 = np.where(lattice == 11)
    x2, y2 = np.where(lattice == 12)
    mol_centers1 = np.array(zip(x1,y1))
    mol_centers2 = np.array(zip(x2,y2))
    mol_centers = np.concatenate((mol_centers1, mol_centers2))
    for coor in mol_centers:
        #print coor
        neigh1 = get_neighbour(coor, lattice.shape[0], 1)
        neigh2 = get_neighbour(coor, lattice.shape[0], 2)
        for i in range(neigh1.shape[0]):
            if lattice[tuple(neigh1[i,:])] == end_group:
                if lattice[tuple(neigh2[i,:])] == 1:
                    count = count + 1
    return count

def num_mbond2(lattice, end_group1, end_group2):
    """ return the number of metal-coordination bond formed by end_group1 and end_group2 in a linear fashion"""
    count = 0
    # get the coordinates of the molecules
    x1, y1 = np.where(lattice == 11)
    x2, y2 = np.where(lattice == 12)
    mol_centers1 = np.array(zip(x1,y1))
    mol_centers2 = np.array(zip(x2,y2))
    mol_centers = np.concatenate((mol_centers1, mol_centers2))
    for coor in mol_centers:
        #print coor
        neigh1 = get_neighbour(coor, lattice.shape[0], 1)
        neigh2 = get_neighbour(coor, lattice.shape[0], 2)
        neigh3 = get_neighbour(coor, lattice.shape[0], 3)
        neigh4 = get_neighbour(coor, lattice.shape[0], 4)
        for i in range(neigh1.shape[0]):
            if lattice[tuple(neigh1[i,:])] == end_group1:
                if lattice[tuple(neigh2[i,:])] == 1:
                    if lattice[tuple(neigh3[i,:])] == end_group2:
                        if lattice[tuple(neigh4[i,:])] in m_center:
                            count = count + 1
    return count


if __name__ == "__main__":
    lattice = np.loadtxt('test.txt',delimiter=',',skiprows=1)
    num_bond = num_mbond(lattice, 3)
    num_bond2 = num_mbond(lattice, 2)
    num_bond3 = num_mbond2(lattice, 3, 2)
    print num_bond, num_bond2, num_bond3



