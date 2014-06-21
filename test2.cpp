// c++ test to load and save txt
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <stdlib.h>
using namespace std;

int main()
{
	int p =0;
	int i = 0;
	int j = 0;
	int a[4][3];
	for(i = 0; i < 4; i = i+1)
	{
		 for(j = 0; j < 3; j = j+1)
		{
			a[i][j] = p;
			p = p + 1;
		}
	}
	
	ofstream file ("lattice.txt");
	string line;
	stringstream linestream;
	for(i = 0; i < 4; i = i+1)
	{
		 for(j = 0; j < 3; j = j+1)
		{
			string line;
			stringstream linestream;
			linestream<<a[i][j];
			linestream>>line;
			cout << a[i][j]<<",";
			if(j !=2)
			{
				file<<line<<",";
				
				}
			else
			{
				file<<line<<"\r\n";
				}
			
		}
	}
	return 0;

}

