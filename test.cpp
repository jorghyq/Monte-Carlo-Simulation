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
	int a[20][20];
	int i = 0;
	int j = 0;
	ifstream file("text.txt");
	string line;
	while(getline(file,line))
	{
		std::stringstream   linestream(line);
		std::string         value;
		while(getline(linestream,value,','))
		{
			//std::cout << "Value(" << value << ")\n";
			stringstream tr;
			tr << value;
			tr >> a[i][j];
			cout << a[i][j] << " " ;
			j = j+1;
			if ((j+1)%20 == 0)
			{
				i = i+1;
				j = 0;
			}
		}
    std::cout << "Line Finished" << std::endl;

}
}
