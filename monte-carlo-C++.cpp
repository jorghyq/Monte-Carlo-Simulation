// C++ version of monte-carlo2.py
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <time.h> 
#include <math.h>
using namespace std;

const int TOTAL_RUN = 100000000;
const int SECOND_LOOP = 1;
const int lattice_size = 50;
const int num_molecule = 72;
const int num_metal = 48;

int gtemp[2];
int lattice[lattice_size][lattice_size];
int molecules[num_molecule][4];
int metals[num_metal][2];
int det_anti(int angle);
int *det_neighbour(int x, int y, int angle);
int is_occupied(int *co, int length);
int is_forbidden(int *co, int length);
int cal_energy_mol(int *ind_molecule);
int cal_energy_metal(int *ind_metal);
void save_to_txt();
void disp_array();
int main()
{
	clock_t start, finish; 	
	srand((unsigned)time(NULL)); 
	memset(lattice, 0, sizeof(lattice[0][0]) * lattice_size * lattice_size);
	memset(molecules, 0, sizeof(molecules[0][0]) * num_molecule * 4);
	memset(metals, 0, sizeof(metals[0][0]) * num_metal * 2);
	for(int i = 0; i < num_molecule; i = i + 1)
	{
		int state = 1;
		while(state == 1)
		{
			int ind_x = rand()%(lattice_size);
			int ind_y = rand()%(lattice_size);
			int angle = rand()%6;
			int angle_anti = det_anti(angle);
			int *points_t1;
		    points_t1 = det_neighbour(ind_x,ind_y,angle);
			int points_temp[4];
			points_temp[0] = ind_x;
			points_temp[1] = ind_y;
			points_temp[2] = angle;
			points_temp[3] = angle_anti;
			if(is_occupied(points_temp,4) == 0)
			{
				
				lattice[points_temp[0]][points_temp[1]] = 3;
				points_t1 = det_neighbour(ind_x,ind_y,angle);
				lattice[points_t1[0]][points_t1[1]] = 2;
				points_t1 = det_neighbour(ind_x,ind_y,angle_anti);
				lattice[points_t1[0]][points_t1[1]] = 2;
				molecules[i][0] = ind_x;
				molecules[i][1] = ind_y;
				molecules[i][2] = angle;
				molecules[i][3] = angle_anti;
				state = 0;
			}
		}
	}
	start = clock();
	for(int i = 0;i<num_metal;i++)
	{
		int state = 1;
		while(state==1)
		{
			int ind_x = rand()%(lattice_size);
			int ind_y = rand()%(lattice_size);
			int points_temp[2] = {ind_x,ind_y};
			if(is_occupied(points_temp,2) == 0 && is_forbidden(points_temp,2) == 0)
			{
				lattice[points_temp[0]][points_temp[1]] = 1;
				metals[i][0] = ind_x;
				metals[i][1] = ind_y;
				state = 0;
				}
			}
		}

	
	int energy_current;
	int energy_new;
	for(int m = 0; m < SECOND_LOOP;m = m + 1)
	{
		for(int l = 0; l < TOTAL_RUN; l = l + 1)
		{
			int ind_ele = rand()%(num_molecule + num_metal); 
			if(ind_ele < num_molecule)
			{
				energy_current = cal_energy_mol(molecules[ind_ele]);
				int *temp;
				lattice[molecules[ind_ele][0]][molecules[ind_ele][1]] = 0;
			    temp = det_neighbour(molecules[ind_ele][0],molecules[ind_ele][1],molecules[ind_ele][2]);
				lattice[temp[0]][temp[1]] = 0;
				temp = det_neighbour(molecules[ind_ele][0],molecules[ind_ele][1],molecules[ind_ele][3]);
				lattice[temp[0]][temp[1]] = 0;
				int state= 1;
				while(state == 1)
				{
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
					int new_angle = rand()%6;
					int new_angle_anti = det_anti(new_angle);
					int new_temp_pos[4] = {new_pos[0],new_pos[1],new_angle,new_angle_anti};
					if ((is_occupied(new_temp_pos,4) == 0) && (is_forbidden(new_temp_pos,4) == 0))
					{
						energy_new = cal_energy_mol(new_temp_pos);
						double p = min(exp(-double(energy_new - energy_current)),double(1));
						if(p > (double)rand()/RAND_MAX)
						{
							lattice[new_pos[0]][new_pos[1]] = 3;
							temp = det_neighbour(new_pos[0],new_pos[1],new_angle);
							lattice[temp[0]][temp[1]] = 2;
							temp = det_neighbour(new_pos[0],new_pos[1],new_angle_anti);
							lattice[temp[0]][temp[1]] = 2;
							molecules[ind_ele][0] = new_pos[0];
							molecules[ind_ele][1] = new_pos[1];
							molecules[ind_ele][2] = new_angle;
							molecules[ind_ele][3] = new_angle_anti;
						}
						else
						{
							lattice[molecules[ind_ele][0]][molecules[ind_ele][1]] = 3;
			        		temp = det_neighbour(molecules[ind_ele][0],molecules[ind_ele][1],molecules[ind_ele][2]);
							lattice[temp[0]][temp[1]] = 2;
							temp = det_neighbour(molecules[ind_ele][0],molecules[ind_ele][1],molecules[ind_ele][3]);
							lattice[temp[0]][temp[1]] = 2;
						}
						state = 0;
					}
				}	
			}
			else
			//if (ind_ele >= num_molecule)
			{
				ind_ele = ind_ele - num_molecule;
				energy_current = cal_energy_metal(metals[ind_ele]);
				lattice[metals[ind_ele][0]][metals[ind_ele][1]] = 0;
				int state = 1;
				while(state == 1)
				{
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
					if(is_occupied(new_pos,2) == 0 && is_forbidden(new_pos,2) == 0)
					{
						energy_new = cal_energy_metal(new_pos);
						double p = min(exp(-double(energy_new - energy_current)),double(1));
						if(p > (double)rand()/RAND_MAX)
						{
							lattice[new_pos[0]][new_pos[1]] = 1;
							metals[ind_ele][0] = new_pos[0];
							metals[ind_ele][1] = new_pos[1];
						}
						else
						{
							lattice[metals[ind_ele][0]][metals[ind_ele][1]] = 1;
						}
						state = 0;
					}
				}
			}
			if((l%(TOTAL_RUN/10)) == 0)
			{
				cout<<"current number: "<< l/(TOTAL_RUN/10)<<endl;
				}
		}
	}
	finish = clock();
	cout<<"costed time: " << (finish-start)/CLOCKS_PER_SEC;
	save_to_txt();
	return 0;
}


int keep_within_n(int n, int var)
{
	int temp;
	if (var > n) {temp = var - n - 1;}
	else if (var < 0) {temp = var + n + 1;}
	else {temp = var;}
	return temp;
	}


int *det_neighbour(int ind_x, int ind_y, int angle)
{
	int i = ind_x;
	int j = ind_y;
	if (angle == 0) {gtemp[0] = i; gtemp[1] = j+1;}
	else if (angle == 1) {gtemp[0] = i-1; gtemp[1] = j+1;}
	else if (angle == 2) {gtemp[0] = i-1; gtemp[1] = j;}
	else if (angle == 3) {gtemp[0] = i; gtemp[1] = j-1;}
	else if (angle == 4) {gtemp[0] = i+1; gtemp[1] = j-1;}
	else {gtemp[0] = i+1; gtemp[1] = j;}
	gtemp[0] = keep_within_n(lattice_size-1, gtemp[0]);
	gtemp[1] = keep_within_n(lattice_size-1, gtemp[1]);
	return gtemp;
}

int det_anti(int angle)
{
	int angle_anti;
	if (angle > 2) angle_anti = angle - 3;
	else angle_anti = angle + 3;
	return angle_anti;
	}
	
int is_occupied(int *co, int length)
{
	if(length == 2)
	{
		if(lattice[co[0]][co[1]] != 0) return 1;
	}
	else
	{
		int *temp;
		if(lattice[co[0]][co[1]] != 0) return 1;
		temp = det_neighbour(co[0],co[1],co[2]);
		if(lattice[temp[0]][temp[1]] != 0) return 1;
		temp = det_neighbour(co[0],co[1],co[3]);
		if(lattice[temp[0]][temp[1]] != 0) return 1;
	}
	return 0;
}

int cal_energy_mol(int *ind_molecule)
{
	int energy = 0;
	int *temp1;
	temp1 = det_neighbour(ind_molecule[0],ind_molecule[1],ind_molecule[2]);
	temp1 = det_neighbour(temp1[0],temp1[1],ind_molecule[2]);
	if(lattice[temp1[0]][temp1[1]] == 1) energy = energy - 10;
	int *temp2;
	temp2 = det_neighbour(ind_molecule[0],ind_molecule[1],ind_molecule[3]);
	temp2 = det_neighbour(temp2[0],temp2[1],ind_molecule[3]);
	if(lattice[temp2[0]][temp2[1]] == 1) energy = energy - 10;
	return energy;
}

int cal_energy_metal(int *ind_metal)
{
	int energy = 0;
	int *temp1;
	for(int i = 0; i < 6; i = i+1)
	{
		temp1 = det_neighbour(ind_metal[0],ind_metal[1],i);
		if(lattice[temp1[0]][temp1[1]] == 2)
		{
			temp1 = det_neighbour(temp1[0],temp1[1],i);
			if(lattice[temp1[0]][temp1[1]] == 3)
			{
				energy = energy - 10;
			}
		}
	}
	return energy;
}

int is_forbidden(int *co, int length)
{
	if(length == 2)
	{
		int count = 0;
		int ind_list[6] = {0,0,0,0,0,0};
		for(int i = 0; i < 6; i++)
		{
			int *temp1;
			temp1 = det_neighbour(co[0],co[1],i);
			if(lattice[temp1[0]][temp1[1]] == 2)
			{
				temp1 = det_neighbour(temp1[0],temp1[1],i);
				if(lattice[temp1[0]][temp1[1]] == 3)
				{
					ind_list[count] = i;
					count = count + 1;
				}
			}
		}
		if(count <= 1) return 0;
		else if(count > 3) return 1;
		else if(count == 2)
		{
			if ((abs(ind_list[0]-ind_list[1]) == 1) || (abs(ind_list[0]-ind_list[1]) == 5)) return 1;
			else return 0;
		}
		else 
		{
			if (((ind_list[2] - ind_list[1]) == 2) and ((ind_list[1] - ind_list[0]) == 2)) return 0;
			else return 1;
		}
	}
	else
	{
		int *temp;
		int temp1[2];
		temp = det_neighbour(co[0],co[1],co[2]);
		temp = det_neighbour(temp[0],temp[1],co[2]);
		temp1[0] = temp[0];
		temp1[1] = temp[1];
		if (lattice[temp[0]][temp[1]] == 1)
		{
			int anti_angle = det_anti(co[2]);
			int anti_angle_nei1 = keep_within_n(5,anti_angle-1);
			temp = det_neighbour(temp[0],temp[1],anti_angle_nei1);
			if (lattice[temp[0]][temp[1]] == 2)
			{
				temp = det_neighbour(temp[0],temp[1],anti_angle_nei1);
				if (lattice[temp[0]][temp[1]] == 3) return 1;
			}
			int anti_angle_nei2 = keep_within_n(5,anti_angle+1);
			temp = det_neighbour(temp1[0],temp1[1],anti_angle_nei2);
			if (lattice[temp[0]][temp[1]] == 2)
			{
				temp = det_neighbour(temp[0],temp[1],anti_angle_nei2);
				if (lattice[temp[0]][temp[1]] == 3) return 1;
			}
		}
		temp = det_neighbour(co[0],co[1],co[3]);
		temp = det_neighbour(temp[0],temp[1],co[3]);
		temp1[0] = temp[0];
		temp1[1] = temp[1];
		if (lattice[temp[0]][temp[1]] == 1)
		{
			int anti_angle = det_anti(co[3]);
			int anti_angle_nei1 = keep_within_n(5,anti_angle-1);
			temp = det_neighbour(temp[0],temp[1],anti_angle_nei1);
			if (lattice[temp[0]][temp[1]] == 2)
			{
				temp = det_neighbour(temp[0],temp[1],anti_angle_nei1);
				if (lattice[temp[0]][temp[1]] == 3) return 1;
			}
			int anti_angle_nei2 = keep_within_n(5,anti_angle+1);
			temp = det_neighbour(temp1[0],temp1[1],anti_angle_nei2);
			if (lattice[temp[0]][temp[1]] == 2)
			{
				temp = det_neighbour(temp[0],temp[1],anti_angle_nei2);
				if (lattice[temp[0]][temp[1]] == 3) return 1;
			}
		}
		return 0;
	}
	
}

void disp_array()
{
	for(int i = 0;i<lattice_size;i++)
	{
		for(int j = 0;j<lattice_size;j++)
		{
			cout<<lattice[i][j]<<" ";
		}
		cout<<endl;
	}
}

void save_to_txt()
{
	ofstream file ("lat.txt");
	string line;
	stringstream linestream;
	for(int i = 0; i < lattice_size; i = i+1)
	{
		 for(int j = 0; j < lattice_size; j = j+1)
		{
			string line;
			stringstream linestream;
			linestream<<lattice[i][j];
			linestream>>line;
			if(j == (lattice_size-1))
			{
				file<<line;
				//cout<<line<<",";
				
				}
			else
			{
				file<<line<<",";
				//cout<<line<<",";
			}
				
			
		}
		file<<"\r\n";
	}
}
