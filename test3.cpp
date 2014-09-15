#include <fstream>
#include <sstream>
#include <iostream>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
//extern char *optarg;
//extern int optind, opterr, optopt;
using namespace std;

const int lattice_size = 10;
int temp[5][2];

int num_molecule = 4;
int num_metal = 4;
int num_total = num_molecule + num_metal;
const int element_num = 150;
int cenergy = 10;
double venergy = 0.5;
int mcenergy = 10;


int lattice[lattice_size][lattice_size];
int lattice_num[lattice_size][lattice_size];
int elements[element_num][2];


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

int direct[5][2] = {{-1,0},{0,1},{1,0},{0,-1},{0,0}};
int ind[5] = {0};


int main(int argc, char *argv[])
{
	ind[0] = 0;
	ind[1] = 1;
	ind[2] = 2;
	ind[3] = 3;
	ind[4] = 4;
	memset(lattice, 0, sizeof(lattice[0][0]) * lattice_size * lattice_size);
	memset(lattice_num, 0, sizeof(lattice_num[0][0]) * lattice_size * lattice_size);
	int ind_x = 2;
	int ind_y = 3;
	int ind_x1 = 0;
	int ind_y1 = 5;
	int ind_mx = 2;
	int ind_my = 5;
	double energy_mol;
	double energy_metal;
	int (*points_t1)[2];
	int (*points_t2)[2];
	int (*points_t3)[2];
	int points_mol[5][2];
	int points_molb[5][2];
	int points_metal[5][2];
	/************************************************************************/
	points_t1 = det_neighbour(ind_x,ind_y,ind,5);	
	p2a(points_t1,points_mol,5);
	set_element(&points_mol[0],5,2,1);
	
	points_t3 = det_neighbour(ind_x1,ind_y1,ind,5);	
	p2a(points_t3,points_molb,5);
	set_element(&points_molb[0],5,2,1);
	
	points_t2 = det_neighbour(ind_mx,ind_my,ind,5);
	p2a(points_t2,points_metal,5);	
	set_element(&points_metal[4],1,1,4);
	
	energy_mol = cal_energy_mol(&points_mol[0],4);
	energy_metal = cal_energy_metal(&points_metal[0],4);
	cout<<" molecule energy: "<<energy_mol<<" metal energy: "<<energy_metal<<endl;
	/*************************************************************************/
	int test1 = is_forbidden(&points_mol[0],4);
	int test2 = is_forbidden(&points_metal[0],1);
	cout<<" molecule is: "<<test1<<" metal is: "<<test2<<endl;
	disp_array(2);
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
		int pos_around2[4][2];
		for (int i=0;i<4;i++)
		{
			pos_around[i][0] = kwn(lattice_size, *co[i] + direct[i][0]);
			pos_around[i][1] = kwn(lattice_size, *(co[i]+1) + direct[i][1]);
			pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
			pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
			if ((lattice_num[pos_around[i][0]][pos_around[i][1]] != 0) && (lattice_num[pos_around[i][0]][pos_around[i][1]] == lattice_num[pos_around2[i][0]][pos_around2[i][1]])) return 1;
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
				//if ((lattice_num[plus1[0]][plus1[1]] == lattice_num[plus2[0]][plus2[1]]) || (lattice_num[minus1[0]][minus1[1]] == lattice_num[minus2[0]][minus2[1]]))
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

//void save_to_txt()
//{
	//string filename;
	//stringstream ss;
	////ss<<total_run<<"-"<<lattice_size<<"-"<<num_molecule<<"-"<<num_metal<<"-"<<cenergy<<"-"<<venergy<<"-"<<mcenergy<<".txt";
	//ss << "results5\\";
	//ss.precision(1);
	//ss.setf(ios::scientific);
	//ss << double(total_run) << "-" << lattice_size << "-" << num_molecule << "-" << num_metal << "-";
	//ss << cenergy << "-" << venergy << "-" << mcenergy << ".txt";
	////ofstream file ("%e-%d-%d-%d-%d-%d-%d.txt", total_run,lattice_size,num_molecule,num_metal,cenergy,venergy,mcenergy);
	//filename = ss.str();
	//cout<<"output to file: "<<filename<<endl;
	//ofstream file(filename.c_str());
	////ofstream file("latt.txt");
	//string line;
	//stringstream linestream;
	//for(int i = 0; i < lattice_size; i = i+1)
	//{
		 //for(int j = 0; j < lattice_size; j = j+1)
		//{
			//string line;
			//stringstream linestream;
			//linestream<<lattice[i][j];
			//linestream>>line;
			//if(j == (lattice_size-1))
			//{
				//file<<line;
				////cout<<line<<",";
				
				//}
			//else
			//{
				//file<<line<<",";
				////cout<<line<<",";
			//}
				
			
		//}
		//file<<"\r\n";
	//}
//}

