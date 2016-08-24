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
const int lattice_size = 80;
const int element_num = 1000;//12 * lattice_size;
int ffn = 1; // number of results filefolder
int num_molecule = 40;
int num_molecule1 = 20;
int num_molecule2 = 20;
int num_metal = 0;
int num_total = num_molecule + num_metal;
double cenergy = 10;
double venergy = 3;
double mcenergy = 10;
int ctemp[3];
int temp[5][3];
int ind[5] = {0,1,2,3,4};
int lattice[lattice_size][lattice_size];
int lattice_num[lattice_size][lattice_size];
int elements[element_num][3];// third column for the angle, 1 means horizontal and 2 means vertical
int direct[5][2] = {{-1,0},{0,1},{1,0},{0,-1},{0,0}};
double energy_sys = 0;
double last_energy =  0;
double p = 0;
double p_temp = 0;
int kwn(int n, int var);
void p2a(int (*ptr)[3], int arr[][3], int length);
int (*det_neighbour(int ind_x, int ind_y, int *ind, int *conf,int length))[3];
int is_occupied(int (*co)[3], int length);
int is_forbidden(int (*co)[3], int length);
void set_element(int (*co)[3], int length, int angle,int ind_ele);
double cal_energy_mol(int (*co)[3], int length);
double cal_energy_metal(int (*co)[3], int length);
void print_array(int (*ar)[3], int length, string name);
void save_to_txt();
void disp_array(int i);
double cal_energy_sys(void);
int *cal_bond_num(void);
/*********************************************************** 
 * molecular macrocycle
 * 9: mol1, light blue
 * 10: mol2, dark blue
 *
 * metal
 * 1: metal
 *
 * molecular components at meso positions
 * 2: 2-fold allowed, red
 * 3: 4-fold allowed, original
 * 4: non-reactive, next to the metal allowed, vdW
 * 5: non-reactive, next to the metal allowed, no vdW
 * 6: non-reactive, next to the metal forbidden, vdW
 * 7: non-reactive, next to the metal forbidden, no vdW\
 * 8: saved
 * **********************************************************/
int mol_conf1[2][5] = {{3,3,3,3,9},{3,3,3,3,9}};
int mol_conf2[2][5] = {{2,2,2,2,10},{2,2,2,2,10}};
int metal_conf[5] = {1,1,1,1,1};


int main(int argc, char *argv[])
{
	// SET THE COMMANDLINE ARGUMENTS
	// a: total_run
	// b: num_mol1
    // c: num_mol2
	// d: num_metal
	// e: cenergy
	// f: venergy
	// g: mcenergy
	int opt;
	const char *optstring = "a:b:c:d:e:f:g:h:";
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
				num_molecule1 = atoi(optarg);
				break;
            case 'c':
                num_molecule2 = atoi(optarg);
			case 'd':
				num_metal = atoi(optarg);
				break;
			case 'e':
				cenergy = atof(optarg);
				break;
			case 'f':
				venergy = atof(optarg);
				break;
			case 'g':
				mcenergy = atof(optarg);
				break;
			case 'h':
				ffn = atoi(optarg);
				break;
			default:
				break;
			}
	}
    num_molecule = num_molecule1 + num_molecule2;
	num_total = num_molecule + num_metal;
	cout<<"Program is initialized with: "<<endl;
	cout<<"total_run = "<< total_run<<", latt_length = "<<lattice_size<<endl;
	cout<<"num mol1 = "<<num_molecule1<< " ,num mol2 = "<<num_molecule2<<", num mol = "<<num_molecule<<", num metals = "<<num_metal<<endl;
	cout<<"cenergy = "<<cenergy<<", venergy = "<<venergy<<", mcenergy = "<<mcenergy<<endl;
	clock_t start, finish; 	
	srand((unsigned)time(NULL)); 
	// INITIALIZE THE LATTICE
	memset(lattice, 0, sizeof(lattice[0][0]) * lattice_size * lattice_size);
	memset(lattice_num, 0, sizeof(lattice_num[0][0]) * lattice_size * lattice_size);
	memset(elements, 0, sizeof(elements[0][0]) * num_total * 3);
	// Distribute the molecules and metals on the lattice
	cout << "Begin to distribute elemenst..." << endl;
    cout<<"molecules1 "<<num_molecule1 << " molecules2 "<<num_molecule2<<endl; 
	for(int i = 0; i < num_molecule + num_metal; i++)
	{
		int state = 1;
		while (state == 1)
		{
			int ind_x = rand()%(lattice_size);
			int ind_y = rand()%(lattice_size);
            int angle = rand()%2+1;// angle = 1 or 2
			//if (angle != 1 and angle != 2) cout<<"ERROR wrong angle number!"<<endl;
            int (*points_t1)[3];
			int points_temp[5][3];
			if (i < num_molecule1)
			{
				points_t1 = det_neighbour(ind_x,ind_y,ind,mol_conf1[angle-1],5);
				p2a(points_t1,&points_temp[0],5);
				//print_array(points_t1,5);
				if ((is_occupied(&points_temp[0],5)) == 0 && (is_forbidden(&points_temp[0],4)) == 0)
				{
					set_element(&points_temp[0],5,angle,i);
                    /* mol1: 5
                     *      565
                     *       5
                     * */
					state = 0;
					//reg2 = reg2 + 1;
				}
				//disp_array(1);
			}
			else if (i < num_molecule)
			{
				points_t1 = det_neighbour(ind_x,ind_y,ind,mol_conf2[angle-1],5);
				p2a(points_t1,&points_temp[0],5);
				//print_array(points_t1,5);
				if ((is_occupied(&points_temp[0],5)) == 0 && (is_forbidden(&points_temp[0],4)) == 0)
				{
					set_element(&points_temp[0],5,angle,i);
                    /*mol2:           3      or           2
                     *      amgle=1  242       angle=2   343
                     *                3                   2
                     * */
					state = 0;
					//reg2 = reg2 + 1;
				}
				//disp_array(1);
			}
			else
			{
				points_t1 = det_neighbour(ind_x,ind_y,ind,metal_conf,5);
				p2a(points_t1,&points_temp[0],5);
			    //print_array(points_t1,4,"points_t1");
				if (is_occupied(&points_temp[4],1) == 0 && is_forbidden(&points_temp[0],1) == 0)
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
    int old_angle;
    int new_angle;
	for(int m = 0; m < SECOND_LOOP;m = m + 1)
	{
		for(int l = 0; l < total_run; l = l + 1)
		{
			int ind_ele = rand()%(num_molecule + num_metal);
            new_angle = rand()%2 + 1;
            //if (new_angle != 1 && new_angle != 2 ) cout<<"ERROR: WRONG angle number!"<<endl;
			//cout << "number of elements " << ind_ele << endl;
			int (*points_t1)[3];
			int points_old[5][3];
			int (*points_t2)[3];
			int points_new[5][3];
            old_angle = elements[ind_ele][2];
			if(ind_ele < num_molecule1)
			{
				//disp_array(1);
                //if (old_angle != 1 && old_angle != 2) cout<<"ERROR: WRONG angle number!"<<endl;
                //cout<<"old angle "<<old_angle<<endl;
                //cout << "number of elements " << ind_ele << endl;
                //cout<<"test0"<<endl;
				points_t1 = det_neighbour(elements[ind_ele][0],elements[ind_ele][1],ind,mol_conf1[old_angle-1],5);
				p2a(points_t1,points_old,5);
				//print_array(points_t1,5,"points_t1");
				//print_array(&points_old[0],5, "points_old");
				int state= 1;
				while (state == 1)
				{
					//cout << "Enters the loop..."<<endl;
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
                    //cout << new_pos[0] <<" "<<new_pos[1]<<endl;
					points_t2 = det_neighbour(new_pos[0],new_pos[1],ind,mol_conf1[new_angle-1],5);
					p2a(points_t2,points_new,5);
					//print_array(points_t2,5);
					if ((is_occupied(&points_new[0],5)) == 0 && (is_forbidden(&points_new[0],4)) == 0)
					{
                        //cout<<"test1"<<endl;
						energy_old = cal_energy_mol(&points_old[0],4);
						energy_new = cal_energy_mol(&points_new[0],4);
                        //cout<<"test2"<<endl;
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<endl;
						p = min(exp(-double(energy_new - energy_old)),double(1));
						p_temp = (rand()+1)/(double(RAND_MAX)+1);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<" probability: "<<p<<" random: "<<temp<<endl;
						//if (p > (double)rand()/RAND_MAX)
						if (p > p_temp)
						{
							//cout<<"energy lowered by "<<energy_new-energy_old<<endl;
							//cout<< "NEW MOLECULE POSITION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
							set_element(&points_old[0],5,0,ind_ele);
						    //print_array(points_t1,5,"points_t1");
							//print_array(&points_old[0],5,"points_old");
							set_element(&points_new[0],5,new_angle,ind_ele);
						}
						state = 0;
						//cout <<ind_ele<< " is done"<<endl;
					}
				}	
			}
			else if(ind_ele < num_molecule)
			{
				//disp_array(1);
                //cout << "number of elements " << ind_ele << endl;
				points_t1 = det_neighbour(elements[ind_ele][0],elements[ind_ele][1],ind,mol_conf2[old_angle-1],5);
				p2a(points_t1,points_old,5);
				//print_array(points_t1,5,"points_t1");
				//print_array(&points_old[0],5, "points_old");
				int state= 1;
				while (state == 1)
				{
					//cout << "Enters the loop..."<<endl;
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
					//int new_angle = rand()%2 + 1;
                    //cout << new_pos[0] <<" "<<new_pos[1]<<endl;
                    points_t2 = det_neighbour(new_pos[0],new_pos[1],ind,mol_conf2[new_angle-1],5);
					p2a(points_t2,points_new,5);
                    //print_array(points_t2,5);
					if ((is_occupied(&points_new[0],5)) == 0 && (is_forbidden(&points_new[0],4)) == 0)
					{
						energy_old = cal_energy_mol(&points_old[0],4);
						energy_new = cal_energy_mol(&points_new[0],4);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<endl;
						p = min(exp(-double(energy_new - energy_old)),double(1));
						p_temp = (rand()+1)/(double(RAND_MAX)+1);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<" probability: "<<p<<" random: "<<temp<<endl;
						//if (p > (double)rand()/RAND_MAX)
						if (p > p_temp)
						{
							//cout<<"energy lowered by "<<energy_new-energy_old<<endl;
							//cout<< "NEW MOLECULE POSITION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
							set_element(&points_old[0],5,0,ind_ele);
							//print_array(points_t1,5,"points_t1");
							//print_array(&points_old[0],5,"points_old");
							set_element(&points_new[0],5,new_angle,ind_ele);
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
				points_t1 = det_neighbour(elements[ind_ele][0],elements[ind_ele][1],ind,metal_conf,5);
				p2a(points_t1,&points_old[0],5);
				while(state == 1)
				{
					int new_pos[2] = {rand()%lattice_size,rand()%lattice_size};
					points_t2 = det_neighbour(new_pos[0],new_pos[1],ind,metal_conf,5);
					p2a(points_t2,&points_new[0],5);
					if (is_occupied(&points_new[4],1) == 0 && is_forbidden(&points_new[0],1) == 0)
					{
						energy_old = cal_energy_metal(&points_old[0],4);
						energy_new = cal_energy_metal(&points_new[0],4);
						//cout<<" old energy: "<<energy_old<<" new energy: "<<energy_new<<endl;
						p = min(exp(-double(energy_new - energy_old)),double(1));
						//cout<<"probability: "<<p<<endl;
						p_temp = (rand()+1)/(double(RAND_MAX)+1);
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

void p2a(int (*ptr)[3], int arr[][3], int length)
{
	for (int i=0;i<length;i++)
	{
		arr[i][0] = *ptr[i];
		arr[i][1] = *(ptr[i]+1);
        arr[i][2] = *(ptr[i]+2);
	}
}
void print_array(int (*ar)[3], int length, string name)
{
	cout<<"Print array: "<<name<<endl;
	for (int i=0;i<length;i++)
	{
		cout<<*ar[i]<<" "<<*(ar[i]+1)<<" "<<*(ar[i]+2)<<endl;
	}

}

// ind is a pointer to a 1d array, where the indexex of wanted points are stored
int (*det_neighbour(int ind_x, int ind_y, int *ind, int *conf, int length))[3]
{
	int ind_num;	
	//cout<<"temp value:"<<endl;
	for (int i=0;i<length;i++)
	{
		ind_num = ind[i];
		temp[i][0] = kwn(lattice_size, ind_x + direct[ind_num][0]);
		temp[i][1] = kwn(lattice_size, ind_y + direct[ind_num][1]);
        temp[i][2] = conf[ind_num];
		//cout<<temp[i][0]<<" "<<temp[i][1]<<" "<<temp[i][2]<<endl;
		}
	return temp;
}

void set_element(int (*co)[3], int length, int angle,int ind_ele)
{
    //int angle = 0;
    if (angle == 0)
    {
        for (int i=0;i<length;i++)
        {
            lattice[*co[i]][*(co[i]+1)] = 0;
        }
    }
    else
    {
        for (int i=0;i<length;i++)
        {
            lattice[*co[i]][*(co[i]+1)] = *(co[i]+2);
            //cout<<"point in "<<*co[i]<<" "<<*(co[i]+1)<<" is set to "<<op<<endl;
        }
    }
    if (length == 1)
    {
        elements[ind_ele][0] = *co[0];
        elements[ind_ele][1] = *(co[0]+1);
        elements[ind_ele][2] = angle;
    }
    else
    {
        elements[ind_ele][0] = *co[length-1];
        elements[ind_ele][1] = *(co[length-1]+1);
        elements[ind_ele][2] = angle;
    }
}


int is_occupied(int (*co)[3], int length)
{
	for (int i=0;i<length;i++)
	{
		if (lattice[*co[i]][*(co[i]+1)] != 0)
		{
			return 1;
		}
	}	
	return 0;
}

double cal_energy_mol(int (*co)[3], int length)
{
	double energy = 0;
	int pos_around[4][3];
	int pos_around2[4][3];
	for (int i=0;i<length;i++)
	{
        //if (lattice[*co[i]][*(co[i]+1)] == 2) continue;
		pos_around[i][0] = kwn(lattice_size, *co[i] + direct[i][0]);
		pos_around[i][1] = kwn(lattice_size, *(co[i]+1) + direct[i][1]);
		pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
		pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
        // coordination bond formed
		if (lattice[pos_around[i][0]][pos_around[i][1]] == 1 && (*(co[i]+2) == 2 || *(co[i]+2) == 3))
		{
			energy = energy - double(cenergy);
		}
        // vdW interaction, vdW reactive endgroup: 2,3,4,6
        if (*(co[i]+2) < 5 || *(co[i]+2) == 6)
        {
            //cout<<"vdW detected"<<endl;
            if ((lattice[pos_around[i][0]][pos_around[i][1]] == 2 || lattice[pos_around[i][0]][pos_around[i][1]] == 3 || lattice[pos_around[i][0]][pos_around[i][1]] == 4 || lattice[pos_around[i][0]][pos_around[i][1]] == 6) && (lattice[pos_around2[i][0]][pos_around2[i][1]] != 9 && lattice[pos_around2[i][0]][pos_around2[i][1]] != 10))
            {
                energy = energy - double(venergy);
                //cout<<"vdW confirmed"<<endl;
            }
        }
	}
	return energy;
}

double cal_energy_metal(int (*co)[3], int length)
{
	double energy = 0;
	int pos_around[4][3];
	int pos_around2[4][3];
	for (int i=0;i<length;i++)
	{
		pos_around[i][0] = *co[i];
		pos_around[i][1] = *(co[i]+1);
		pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
		pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
		if (lattice[pos_around[i][0]][pos_around[i][1]] == 2 || lattice[pos_around[i][0]][pos_around[i][1]] == 3)
		{
			if (lattice[pos_around2[i][0]][pos_around2[i][1]] > 8)
			{
				energy = energy - double(mcenergy);
                //cout<<"coordination detected!"<<endl;
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
	int (*points_t1)[3];
    int angle = 0;
	for (int i=0;i<num_total;i++)
	{
		energy_temp = 0;
        angle = elements[i][2];
		if (i < num_molecule1)
		{
            points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,mol_conf1[angle-1],5);
			energy_temp = cal_energy_mol(points_t1,4);
		}
        else if(i < num_molecule)
        {
            points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,mol_conf2[angle-1],5);
			energy_temp = cal_energy_mol(points_t1,4);
        }
		else
		{
			points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,metal_conf,5);
			energy_temp = cal_energy_metal(points_t1,4);
		}
		energy = energy + energy_temp;
		//cout<<"energy change: "<<energy_temp<<" energy total: "<<energy<<endl;
	}
	return energy/2.0;	
}

int is_forbidden(int (*co)[3], int length)
{
	
	if(length == 1)
	{
		int pos_around[4][3];
		int pos_around2[4][3];
		int count[4] = {0,0,0,0};
		int count_num = 0;
		//pos_around = det_neighbour(*co[0],*(co[0]+1),ind,4);
		for (int i=0;i<4;i++)
		{
			pos_around[i][0] = *co[i];
			pos_around[i][1] = *(co[i]+1);
			pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
			pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
			// metals can not form cluster
            if (lattice[pos_around[i][0]][pos_around[i][1]] == 1) return 1;
            // endgroup 6 and 7 are forbidden for coordination
            else if (lattice[pos_around[i][0]][pos_around[i][1]] == 6 || lattice[pos_around[i][0]][pos_around[i][1]] == 7)
            {
                if (lattice[pos_around2[i][0]][pos_around2[i][1]] > 8)  return 1;
            }
            // check if the linear 2-fold coordination condition is violated
			else if (lattice[pos_around[i][0]][pos_around[i][1]] == 2)
			{
				if (lattice[pos_around2[i][0]][pos_around2[i][1]] > 8)
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
            //cout<<"Check if forbidden"<<endl;
			if ((count[0] == 1) && (count[2] == 1)) return 0;
			else if ((count[1] == 1) && (count[3] == 1)) return 0;
			else return 1;
		}
		else return 0;
	}
	else
	{
		int pos_around[4][3];
		for (int i=0;i<4;i++)
		{
            // if the meso substituent is only allowed for two-fold coordination
            pos_around[i][0] = kwn(lattice_size, *co[i] + direct[i][0]);
            pos_around[i][1] = kwn(lattice_size, *(co[i]+1) + direct[i][1]);
            if (*(co[i]+2) == 2)
            {
                //pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
                //pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
                int plus1[2] = {0};
                int plus2[2] = {0};
                int minus1[2] = {0};
                int minus2[2] = {0};
				if (lattice[pos_around[i][0]][pos_around[i][1]] == 1)
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
                    // 4-fold not allowed, here 2-fold coordination could also
                    // have 4-fold at its neighbour coordination positions
                    if (lattice[plus1[0]][plus1[1]] == 2 && lattice[plus2[0]][plus2[1]] > 8)
                    {
                        return 1;
                    }
                    if (lattice[minus1[0]][minus1[1]] == 2 && lattice[minus2[0]][minus2[1]] > 8)
                    {
                        return 1;
                    }
                }
            }
            // if the meso substituent is non active, just continue
            else if (*(co[i]+2) == 6 || *(co[i]+2) == 7)
            {
				if (lattice[pos_around[i][0]][pos_around[i][1]] == 1)
                {
                    return 1;
                }
            }
            //else if (*(co[i]+2) == 7) continue;
            //else if (*(co[i]+2) == 5) continue;
            // if 4-fold no restrictions, *(co[i]+2) == 5
            else continue;   
		}
		return 0;
	}
}

int *cal_bond_num(void)
{
	int cbond = 0;
	int vbond = 0;
	int (*points_t1)[3];
	int pos_around[4][3];
	int pos_around2[4][3];
    int angle = 0;
	for (int i=0;i<num_molecule;i++)
	{
        angle = elements[i][2];
        if (i < num_molecule1)
		{
            points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,mol_conf1[angle-1],5);
		}
        else if(i < num_molecule)
        {
            points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,mol_conf2[angle-1],5);
        }
		points_t1 = det_neighbour(elements[i][0],elements[i][1],ind,mol_conf2[angle-1],5);
		for (int i=0;i<4;i++)
		{
			pos_around[i][0] = kwn(lattice_size,*points_t1[i] + direct[i][0]);
			pos_around[i][1] = kwn(lattice_size,*(points_t1[i]+1) + direct[i][1]);
			pos_around2[i][0] = kwn(lattice_size, pos_around[i][0] + direct[i][0]);
			pos_around2[i][1] = kwn(lattice_size, pos_around[i][1] + direct[i][1]);
			if (lattice[pos_around[i][0]][pos_around[i][1]] == 1 && *(points_t1[i]+2) < 4)
			{
				cbond = cbond + 1;
			}
			else if (lattice[pos_around[i][0]][pos_around[i][1]] != 0)
			{
				if (lattice[pos_around[i][0]][pos_around[i][1]] != lattice[pos_around2[i][0]][pos_around2[i][1]])
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
	//ss << "D:\\Dropbox\\Project\\python\\Monte-Carlo-Simulation\\results";
	ss << "/home/jorghyq/Dropbox/Project/python/Monte-Carlo-Simulation/results";
	ss << ffn << "/";
	ss.precision(1);
	ss.setf(ios::scientific);
	ss << double(total_run) << "_" << lattice_size << "_" << num_molecule1 << "_" << num_molecule2 << "_";
	ss << num_metal << "_"<< cenergy << "_" << venergy << "_" << mcenergy << ".txt";
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
	/*for(int i = 0; i < num_total; i = i+1)
	{
		 for(int j = 0; j < 3; j = j+1)
		{
			
			//string line;
			//stringstream linestream;
			linestream<<elements[i][j];
			linestream>>line;
			linestream.clear();
			if(j == 2)
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
	}*/
    /*for(int i = 0; i < lattice_size; i = i+1)
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
	}*/
}
