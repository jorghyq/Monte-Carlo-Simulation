// C++ version of monte-carlo2.py
// Without using lattice
// 
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <time.h> 
#include <math.h>
#include <unistd.h>
using namespace std;

int total_run = 100000000;
const int SECOND_LOOP = 1;
const int lattice_size = 100;
const int element_num = 800;//12 * lattice_size;
int num_molecule = 50;
int num_metal = 50;
int num_total = num_molecule + num_metal;
int cenergy = 10;
// needed Functions 
int kwn(int n, int var);
int (*det_shape(int ind_x, int ind_y, int *ind, int length))[2];
double cal_energy(int (*co)[2], int length);
void set_element(int (*co)[2], int length, int op,int ind_ele);
void print_array(int (*ar)[2], int length, string name);
void save_to_txt();
void disp_array(int i);

// Main functions
void make_move();
double cal_energy(int (*co)[2], int length);

int main()
{
	// Initialize the elements
	memset(elements, 0, sizeof(elements[0][0]) * num_total * 4);
	// Initialize the shape of the elements
	int mlenght = 5;
	int 
	}
