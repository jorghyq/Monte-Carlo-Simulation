// C++ version of monte-carlo2.py
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


int total_run = 100000;
const int SECOND_LOOP = 1;
const int lattice_size = 100;
const int element_num = 1000;//12 * lattice_size;
int num_molecule = 200;
int num_metal = 400;
int num_total = num_molecule + num_metal;
int cenergy = 10;
double venergy = 3;
int mcenergy = 10;
int ctemp[2];
int temp[5][2];
int ind[5] = {0};
int lattice[lattice_size][lattice_size];
int lattice_num[lattice_size][lattice_size];
int elements[element_num][2];
int direct[5][2] = {{-1,0},{0,1},{1,0},{0,-1},{0,0}};
double energy_sys = 0;
double last_energy =  0;
double p = 0;
double p_temp = 0;
int kwn(int n, int var);
void p2a(int (*ptr)[2], int arr[][2], int length);
int (*det_neighbour(int ind_x, int ind_y, int *ind, int length))[2];
int is_occupied(int (*co)[2], int length);
int is_forbidden(int (*co)[2], int length);
void set_element(int (*co)[2], int length, int op,int ind_ele);
double cal_energy_mol(int (*co)[2], int length);
double cal_energy_metal(int (*co)[2], int length);
void print_array(int (*ar)[2], int length, string name);
void save_to_txt();
void disp_array(int i);
double cal_energy_sys(void);
int *cal_bond_num(void);

int main(int argc, char *argv[])
{
	// SET THE COMMANDLINE ARGUMENTS
	// a: total_run
	// b: num_mol
	// c: num_metal
	// d: cenergy
	// e: venergy
	// f: mcenergy
	int opt;
	const char *optstring = "a:b:c:d:e:f:";
	while((opt = getopt(argc, argv, optstring)) != -1)
	{
		switch(opt)
		{
			case 'a':
				total_run = atoi(optarg);
				break;
			//case 'b':
			//	lattice_size = atoi(optarg);
			//	break;
			case 'b':
				num_molecule = atoi(optarg);
				break;
			case 'c':
				num_metal = atoi(optarg);
				break;
			case 'd':
				cenergy = atoi(optarg);
				break;
			case 'e':
				venergy = atof(optarg);
				break;
			case 'f':
				mcenergy = atoi(optarg);
				break;	
			default:
				break;
			}
	}
	num_total = num_molecule + num_metal;
	cout<<"Program is initialized with: "<<endl;
	cout<<"total_run = "<< total_run<<", latt_length = "<<lattice_size<<endl;
	cout<<"num of molecule = "<<num_molecule<<", num of metals = "<<num_metal<<endl;
	cout<<"cenergy = "<<cenergy<<", venergy = "<<venergy<<", mcenergy = "<<mcenergy<<endl;
	clock_t start, finish; 	
	srand((unsigned)time(NULL)); 
	// INITIALIZE THE LATTICE
	memset(lattice, 0, sizeof(lattice[0][0]) * lattice_size * lattice_size);
	memset(lattice_num, 0, sizeof(lattice_num[0][0]) * lattice_size * lattice_size);
	memset(elements, 0, sizeof(elements[0][0]) * num_total * 4);
	// Distribute the molecules and metals on the lattice
	cout << "Begin to distribute elemenst..." << endl;
	ind[0] = 0;
	ind[1] = 1;
	ind[2] = 2;
	ind[3] = 3;
	ind[4] = 4;
	
	for(int i = 0; i < num_molecule + num_metal; i++)
	{
		int state = 1;
		while (state == 1)
		{
			int ind_x = rand()%(lattice_size);
			int ind_y = rand()%(lattice_size);
			int (*points_t1)[2];
			int points_temp[5][2];
			if (i < num_molecule)
			{
				points_t1 = det_neighbour(ind_x,ind_y,ind,5);
				p2a(points_t1,&points_temp[0],5);
				//print_array(points_t1,5);
				if ((is_occupied(&points_temp[0],5)) == 0)// && (is_forbidden(&points_temp[0],4)) == 0)
				{
					set_element(&points_temp[0],5,2,i);
					state = 0;
					//reg2 = reg2 + 1;
				}
				//disp_array(1);
			}
			else
			{
				points_t1 = det_neighbour(ind_x,ind_y,ind,5);
				p2a(points_t1,&points_temp[0],5);
				//print_array(points_t1,4);
				if (is_occupied(&points_temp[0],1) == 0)// && is_forbidden(&points_temp[0],1) == 0)
				{
					set_element(&points_temp[4],1,1,i);
					state = 0;
					//reg = reg + 1;
				}
				//disp_array(1);		
			}
		}
		//cout << "The" << i << "is done..."<< endl;
	}
	cout << "Elements are distributed..." << endl;
	start = clock();
	//disp_array(1);
	// Begin to simulatte
	cout << "Begin to simulate..." << endl;
	double energy_old;
	double energy_new;
	for(int m = 0; m < SECOND_LOOP;m = m + 1)
	{
		for(int l = 0; l < total_run; l = l + 1)
		{
			int ind_ele = rand()%(num_molecule + num_metal);
			//cout << "number of elements " << ind_ele << endl;
			int (*points_t1)[2];
			int points_old[5][2];
			int (*points_t2)[2];
			int points_new[5][2];
			if(ind_ele < num_molecule)
			{
				//disp_array(1);
				//cout << "number of elements " << ind_ele << endl;
				points_t1 = det_neighbour(elements[ind_ele][0],elements[ind_ele][1],ind,5);
				p2a(points_t1,points_old,5);
				//print_array(points_t1,5,"points_t1");
				//print_array(&points_old[0],5, "points_old");
				int state= 1;
				while (state == 1)
				{
					//cout << "Enters the loop..."<<endl;
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
					//cout << new_pos[0] <<" "<<new_pos[1]<<endl;
					points_t2 = det_neighbour(new_pos[0],new_pos[1],ind,5);
					p2a(points_t2,points_new,5);
					//print_array(points_t2,5);
					if ((is_occupied(&points_new[0],5)) == 0)// && (is_forbidden(&points_new[0],4)) == 0)
					{
						energy_old = cal_energy_mol(&points_old[0],4);
						energy_new = cal_energy_mol(&points_new[0],4);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<endl;
						p = min(exp(-double(energy_new - energy_old)),double(1));
						p_temp = (rand()+1)/double(RAND_MAX+1);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<" probability: "<<p<<" random: "<<temp<<endl;
						//if (p > (double)rand()/RAND_MAX)
						if (p > p_temp)
						{
							//cout<<"energy lowered by "<<energy_new-energy_old<<endl;
							//cout<< "NEW MOLECULE POSITION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
							set_element(&points_old[0],5,0,ind_ele);
							//print_array(points_t1,5,"points_t1");
							//print_array(&points_old[0],5,"points_old");
							set_element(&points_new[0],5,2,ind_ele);
						}
						state = 0;
						//cout <<ind_ele<< " is done"<<endl;
					}
				}	
			}
			else
			//if (ind_ele > num_molecule-1)
			{
				int state = 1;
				points_t1 = det_neighbour(elements[ind_ele][0],elements[ind_ele][1],ind,5);
				p2a(points_t1,&points_old[0],5);
				while(state == 1)
				{
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
					points_t2 = det_neighbour(new_pos[0],new_pos[1],ind,5);
					p2a(points_t2,&points_new[0],5);
					if (is_occupied(&points_new[0],1) == 0)// && is_forbidden(&points_new[0],1) == 0)
					{
						energy_old = cal_energy_metal(&points_old[0],4);
						energy_new = cal_energy_metal(&points_new[0],4);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<endl;
						p = min(exp(-double(energy_new - energy_old)),double(1));
						//cout<<"probability: "<<p<<endl;
						p_temp = (rand()+1)/double(RAND_MAX+1);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<" probability: "<<p<<" random: "<<temp<<endl;
						//if(p > (double)rand()/RAND_MAX)
						if(p > p_temp)
						{
							//cout<<"energy lowered by "<<energy_new-energy_old<<endl;
							//cout<< "NEW METAL POSITION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
							set_element(&points_old[4],1,0,ind_ele);
							set_element(&points_new[4],1,1,ind_ele);
						}
						state = 0;
						//cout<<ind_ele<<" is done"<<endl;
					}
				}
			}
			//disp_array(1);
			//cout << "current " << l << endl;
			
			energy_sys = cal_energy_sys();
			//if ((last_energy - energy_sys) < -5)
			//{
			//	cout<<" last energy: "<<last_energy<<" new energy: "<<energy_sys<<endl;
			//	cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<" probability: "<<p<<" random: "<<p_temp<<endl;

			//	} 
			last_energy = energy_sys;
			//cout<<" system energy: "<<energy_sys<<endl;
			if((l%(total_run/10)) == 0)
			{
				finish = clock();
				energy_sys = cal_energy_sys();
				cout<<"current number: "<< l/(total_run/10)<<",time: "<<(finish-start)/CLOCKS_PER_SEC<<" system energy: "<<energy_sys<<endl;
			}
		}
	}
	finish = clock();
	energy_sys = cal_energy_sys();
	cout<<"costed time: " << (finish-start)/CLOCKS_PER_SEC<<" system energy: "<<energy_sys<<endl;
	save_to_txt();
	int *bond_num;
	bond_num = cal_bond_num();
	cout<<"cbond: "<<bond_num[0]<<" vbond: "<<bond_num[1]<<endl;
	//disp_array(1);
	//cout<<endl;
	//disp_array(2);
	return 0;
}


int kwn(int n, int var)
{
	int temp;
	if (var > n-1) {temp = var - n;}
	else if (var < 0) {temp = var + n;}
	else {temp = var;}
	return temp;
	}

void p2a(int (*ptr)[2], int arr[][2], int length)
{
	for (int i=0;i<length;i++)
	{
		arr[i][0] = *ptr[i];
		arr[i][1] = *(ptr[i]+1);
	}
}
void print_array(int (*ar)[2], int length, string name)
{
	cout<<"Print array: "<<name<<endl;
	for (int i=0;i<length;i++)
	{
		cout<<*ar[i]<<" "<<*(ar[i]+1)<<endl;
	}

}

// ind is a pointer to a 1d array, where the indexex of wanted points are stored
int (*det_neighbour(int ind_x, int ind_y, int *ind, int length))[2]
{
	int ind_num;	
	//cout<<"temp value:"<<endl;
	for (int i=0;i<length;i++)
	{
		ind_num = ind[i];
		temp[i][0] = kwn(lattice_size, ind_x + direct[ind_num][0]);
		temp[i][1] = kwn(lattice_size, ind_y + direct[ind_num][1]);
		//cout<<temp[i][0]<<" "<<temp[i][1]<<endl;
		}
	return temp;
}

void set_element(int (*co)[2], int length, int op,int ind_ele)
{
	for (int i=0;i<length;i++)
	{
		lattice[*co[i]][*(co[i]+1)] = op;
		if (length == 5 && i == 4 && op != 0) lattice[*co[i]][*(co[i]+1)] = op+1;
		//cout<<"point in "<<*co[i]<<" "<<*(co[i]+1)<<" is set to "<<op<<endl;
		if (op == 0)
		{
			lattice_num[*co[i]][*(co[i]+1)] = op;
			}
		else
		{
			if (ind_ele < num_molecule)
			{
				lattice_num[*co[i]][*(co[i]+1)] = ind_ele+1;
				}
			else
			{
				lattice_num[*co[i]][*(co[i]+1)] = num_total+1;
				}
			}
		}
	if (length == 1)
	{
		elements[ind_ele][0] = *co[0];
		elements[ind_ele][1] = *(co[0]+1);
	}
	else
	{
		elements[ind_ele][0] = *co[length-1];
		elements[ind_ele][1] = *(co[length-1]+1);
	}
}


int is_occupied(int (*co)[2], int length)
{
	if (length == 1)
	{
		if (lattice[*co[4]][*(co[4]+1)] != 0)
			return 1;
	}
	for (int i=0;i<length;i++)
	{
		if (lattice[*co[i]][*(co[i]+1)] != 0)
		{
			return 1;
		}
	}	
	return 0;
}

double cal_energy_mol(int (*co)[2], int length)
{
	double energy = 0;
	int pos_around[4][2];
	int pos_around2[4][2];
	for (int i=0;i<length;i++)
	{
		pos_around[i][0] = kwn(lattice_size, *co[i] + direct[i][0]);
		pos_around[i][1] = kwn(lattice_size, *(co[i]+1) + direct[i][1]);
		pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
		pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
		if (lattice_num[pos_around[i][0]][pos_around[i][1]] == (num_total+1))
		{
			energy = energy - double(cenergy);
			
			}
		else if (lattice_num[pos_around[i][0]][pos_around[i][1]] != 0)
		{
			if (lattice_num[pos_around[i][0]][pos_around[i][1]] != lattice_num[pos_around2[i][0]][pos_around2[i][1]])
			{
				energy = energy - double(venergy);
				}
			}
		}
	return energy;
}

double cal_energy_metal(int (*co)[2], int length)
{
	double energy = 0;
	int pos_around[4][2];
	int pos_around2[4][2];
	for (int i=0;i<length;i++)
	{
		pos_around[i][0] = *co[i];
		pos_around[i][1] = *(co[i]+1);
		pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
		pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
		if ((lattice_num[pos_around[i][0]][pos_around[i][1]] != 0) && (lattice_num[pos_around[i][0]][pos_around[i][1]] != num_total+1))
		{
			if (lattice_num[pos_around[i][0]][pos_around[i][1]] == lattice_num[pos_around2[i][0]][pos_around2[i][1]])
			{
				energy = energy - double(mcenergy);
				}
			}
		}
	return energy;
}

double cal_energy_sys(void)
{
	double energy = 0;
	double energy_temp = 0;
	//int cbond_count = 0;
	//int vbond_count = 0;
	//int temp = 0;
	int (*points_t1)[2];
	for (int i=0;i<num_total;i++)
	{
		energy_temp = 0;
		if (i < num_molecule)
		{
			points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,5);
			energy_temp = cal_energy_mol(points_t1,4);
			//if (energy_temp/cenergy != 0)
			//{
			//	cbond_count = cbond_count+
			//	}
			}
		else
		{
			points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,5);
			energy_temp = cal_energy_metal(points_t1,4);
			
			}
		energy = energy + energy_temp;
		//cout<<"energy change: "<<energy_temp<<" energy total: "<<energy<<endl;
	}
	return energy/2.0;	
}

int is_forbidden(int (*co)[2], int length)
{
	
	if(length == 1)
	{
		int pos_around[4][2];
		int pos_around2[4][2];
		int count[4] = {0};
		int count_num = 0;
		//pos_around = det_neighbour(*co[0],*(co[0]+1),ind,4);
		for (int i=0;i<4;i++)
		{
			pos_around[i][0] = *co[i];
			pos_around[i][1] = *(co[i]+1);
			pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
			pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
			if (lattice_num[pos_around[i][0]][pos_around[i][1]] == (num_total+1)) return 1;
			else if ((lattice_num[pos_around[i][0]][pos_around[i][1]] != 0) && (lattice_num[pos_around[i][0]][pos_around[i][1]] != (num_total+1)))
			{
				if (lattice_num[pos_around[i][0]][pos_around[i][1]] == lattice_num[pos_around2[i][0]][pos_around2[i][1]])
				{
					count[i] = 1;
					count_num = count_num +1;
					}
				}
			}
		//cout<<"num of count: "<<count_num<<" "<<count[0]<<" "<<count[1]<<" "<<count[2]<<" "<<count[3]<<endl;
		if (count_num > 2) return 1;
		else if (count_num == 2)
		{
			if ((count[0] == 1) && (count[2] == 1)) return 0;
			else if ((count[1] == 1) && (count[3] == 1)) return 0;
			else return 1;
			}
		else return 0;
	}
	else
	{
		int pos_around[4][2];
		//int pos_around2[4][2];
		for (int i=0;i<4;i++)
		{
			pos_around[i][0] = kwn(lattice_size, *co[i] + direct[i][0]);
			pos_around[i][1] = kwn(lattice_size, *(co[i]+1) + direct[i][1]);
			//pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
			//pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
			//if ((lattice_num[pos_around[i][0]][pos_around[i][1]] != 0) && (lattice_num[pos_around[i][0]][pos_around[i][1]] == lattice_num[pos_around2[i][0]][pos_around2[i][1]])) return 1;
			int plus1[2] = {0};
			int plus2[2] = {0};
			int minus1[2] = {0};
			int minus2[2] = {0};
			if (lattice_num[pos_around[i][0]][pos_around[i][1]] == (num_total+1))
			{
				int plus = kwn(4,i+1);
				int minus = kwn(4,i-1);
				//cout<<"i: "<<i<<" plus: "<<plus<<" minus: "<<minus<<endl;
				
				//cout<<"minus"<<minus<<endl;
				// points 
				plus1[0] = kwn(lattice_size, pos_around[i][0] + direct[plus][0]);
				plus1[1] = kwn(lattice_size, pos_around[i][1] + direct[plus][1]);
				//cout<<"plus1x: "<<plus1[0]<<" plus1y: "<<plus1[1]<<endl;
				plus2[0] = kwn(lattice_size, plus1[0] + direct[plus][0]);
				plus2[1] = kwn(lattice_size, plus1[1] + direct[plus][1]);
				//cout<<"plus2x: "<<plus2[0]<<" plus2y: "<<plus2[1]<<endl;
				// points
				minus1[0] = kwn(lattice_size, pos_around[i][0] + direct[minus][0]);
				minus1[1] = kwn(lattice_size, pos_around[i][1] + direct[minus][1]);
				minus2[0] = kwn(lattice_size, minus1[0] + direct[minus][0]);
				minus2[1] = kwn(lattice_size, minus1[1] + direct[minus][1]);
				if ((lattice_num[plus1[0]][plus1[1]] != 0) && (lattice_num[plus1[0]][plus1[1]] == lattice_num[plus2[0]][plus2[1]]))
				{
					return 1;
				}
				if ((lattice_num[minus1[0]][minus1[1]] != 0) && (lattice_num[minus1[0]][minus1[1]] == lattice_num[minus2[0]][minus2[1]]))
				{
					return 1;
				}
						
			}
				//else return 0;
		}
		return 0;
	}
	
}

int *cal_bond_num(void)
{
	int cbond = 0;
	int vbond = 0;
	int (*points_t1)[2];
	int pos_around[4][2];
	int pos_around2[4][2];
	for (int i=0;i<num_molecule;i++)
	{
		points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,5);
		for (int i=0;i<4;i++)
		{
			pos_around[i][0] = kwn(lattice_size,*points_t1[i] + direct[i][0]);
			pos_around[i][1] = kwn(lattice_size,*(points_t1[i]+1) + direct[i][1]);
			pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
			pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
			if (lattice_num[pos_around[i][0]][pos_around[i][1]] == (num_total+1))
			{
				cbond = cbond + 1;
				}
			else if (lattice_num[pos_around[i][0]][pos_around[i][1]] != 0)
			{
				if (lattice_num[pos_around[i][0]][pos_around[i][1]] != lattice_num[pos_around2[i][0]][pos_around2[i][1]])
				{
					vbond = vbond + 1;
					}
				}
			}
		}
	ctemp[0] = cbond;
	ctemp[1] = vbond/2;
	return ctemp;
}

void disp_array(int i)
{
	int nmol = 0;
	int nmet = 0;
	if (i == 1)
	{
		for(int i = 0;i<lattice_size;i++)
		{
			for(int j = 0;j<lattice_size;j++)
			{
				if (lattice[i][j] == 1) nmet = nmet + 1;
				else if (lattice[i][j] == 2) nmol = nmol + 1;
				cout<<lattice[i][j]<<" ";
			}
			cout<<endl;
		}
		//cout<<endl;
		cout << "number of metal: " << nmet << " number of molecules: "<< nmol/5 <<endl;
	}
	else
	{
		for(int i = 0;i<lattice_size;i++)
		{
			for(int j = 0;j<lattice_size;j++)
			{
				cout<<lattice_num[i][j]<<" ";
			}
			cout<<endl;
		}
	}
}

void save_to_txt()
{
	string filename;
	stringstream ss;
	//ss<<total_run<<"-"<<lattice_size<<"-"<<num_molecule<<"-"<<num_metal<<"-"<<cenergy<<"-"<<venergy<<"-"<<mcenergy<<".txt";
	ss << "results9\\";
	ss.precision(1);
	ss.setf(ios::scientific);
	ss << double(total_run) << "-" << lattice_size << "-" << num_molecule << "-" << num_metal << "-";
	ss << cenergy << "-" << venergy << "-" << mcenergy << ".txt";
	filename = ss.str();
	cout<<"output to file: "<<filename<<endl;
	ofstream file(filename.c_str());
	//ofstream file("latt.txt");
	string line;
	stringstream linestream;
	int *bond_num;
	double enersys;
	bond_num = cal_bond_num();
	linestream<<bond_num[0];
	linestream>>line;
	file<<line<<",";
	linestream.clear();
	linestream<<bond_num[1];
	linestream>>line;
	file<<line<<",";
	linestream.clear();
	enersys = cal_energy_sys();
	linestream<<enersys;
	linestream>>line;
	file<<line<<",";
	linestream.clear();
	linestream<<lattice_size;
	linestream>>line;
	file<<line<<"\n";
	linestream.clear();
	//cout<<"cbond: "<<bond_num[0]<<" vbond: "<<bond_num[1]<<endl;
	for(int i = 0; i < lattice_size; i = i+1)
	{
		 for(int j = 0; j < lattice_size; j = j+1)
		{
			
			//string line;
			//stringstream linestream;
			linestream<<lattice[i][j];
			linestream>>line;
			linestream.clear();
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
	file<<"\r\n";
	for(int i = 0; i < lattice_size; i = i+1)
	{
		 for(int j = 0; j < lattice_size; j = j+1)
		{
			
			//string line;
			//stringstream linestream;
			linestream<<lattice_num[i][j];
			linestream>>line;
			linestream.clear();
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
