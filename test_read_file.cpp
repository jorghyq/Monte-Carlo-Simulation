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
int output[2][5];
int (*read_conf(int ind))[5];
void read_lattice_from_txt(const char* file);
void read_element_from_txt(const char* file);
const int lattice_size = 8;
void print_array(int (*ar)[3], int length, string name);
void disp_array();
int lattice[lattice_size][lattice_size];
int elements[10][3];

int main ()
{
    memset(lattice, 0, sizeof(lattice[0][0]) * lattice_size * lattice_size);
	//memset(lattice_num, 0, sizeof(lattice_num[0][0]) * lattice_size * lattice_size);
	memset(elements, 0, sizeof(elements[0][0]) * 10 * 3);
    int array[2][5];
    int (*temp)[5];
    temp = read_conf(1);
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<5;j++)
            cout<<*(temp[i]+j)<<"  ";
        cout<<endl;
    }
    string file = "test";
    string file_lattice = file + ".txt";
    string file_element = file + "_element.txt";
    const char * f_lattice = file_lattice.c_str();
    const char * f_element = file_element.c_str();
    cout << f_lattice<<endl;
    cout << f_element<<endl;
    //read_lattice_from_txt(f_lattice);
    read_element_from_txt(f_element);
    //disp_array();
    print_array(elements, 10, "element");
    return 0;
}

int (*read_conf(int ind))[5]
{
    fstream infile;
    int i=0;
    int j=0;
    char cNum[256] ;
    infile.open("mol1.txt", ifstream::in);
    if (infile.is_open())
    {
         while (infile.good())
             {
                  infile.getline(cNum, 256, ',');
                  output[j][i]= atoi(cNum) ;
                  cout<<j<<"  "<<i<<"  "<<output[j][i]<<" ";
                  cout <<"string "<< cNum <<endl;
                  i++;
                  if (i == 5)
                  {
                      j++;
                      i=0;
                  }

             }
             infile.close();
    }
    else
    {
        cout << "Error opening file";
    }
    return output;
}

void read_lattice_from_txt(const char *file)
{
    fstream infile;
    int i=0;
    int j=0;
    char cNum[256];
    int skip = 4;
    infile.open(file, ifstream::in);
    if (infile.is_open())
    {
         infile.getline(cNum, 256);
         cout<<"header "<<cNum<<endl;
         while (infile.good())
             {
                 if(i == lattice_size-1)
                 {
                     infile.getline(cNum, 256);
                 }
                 else
                 {
                     infile.getline(cNum, 256, ',');
                 }
                 lattice[j][i]= atoi(cNum);
                 cout<<j<<"  "<<i<<"  "<<lattice[j][i]<<" ";
                 cout <<"string "<< cNum <<endl;
                 if(i == lattice_size-1)
                 {
                     j++;
                     i=0;
                     cout << "reset"<<endl;
                 }
                 else
                 {
                     i++;
                 }
             }
             infile.close();
    }
    else
    {
        cout << "Error opening file";
    }

}

void read_element_from_txt(const char *file)
{
    fstream infile;
    int i=0;
    int j=0;
    char cNum[256] ;
    infile.open(file, ifstream::in);
    if (infile.is_open())
    {
         while (infile.good())
             {
                  if(i == 2)
                     {
                         infile.getline(cNum, 256);
                     }
                     else
                     {
                         infile.getline(cNum, 256, ',');
                     }
                     elements[j][i]= atoi(cNum);
                     cout<<j<<"  "<<i<<"  "<<elements[j][i]<<" ";
                     cout <<"string "<< cNum <<endl;
                     if(i == 2)
                     {
                         j++;
                         i=0;
                         cout << "reset"<<endl;
                     }
                     else
                     {
                         i++;
                     }
             }
             infile.close();
    }
    else
    {
        cout << "Error opening file";
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
    cout<<endl;
}
