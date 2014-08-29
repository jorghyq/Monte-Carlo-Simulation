#include <fstream>
#include <sstream>
#include <iostream>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
//extern char *optarg;
//extern int optind, opterr, optopt;
using namespace std;
const int num1 = 2;
const int num2 = 3;
int testprint(int (*coor)[2])
{
	cout << *coor[1] << " " << *(coor[1]+1) << endl;
	
	return 0;
}
int main(int argc, char *argv[])
{
	int opt;
	int test[4];
	//int i = 0;
	string temp;
	
	const char *optstring = "a:b:c:d:";
	while((opt = getopt(argc, argv, optstring)) != -1)
	{
		switch(opt)
		{
			case 'a':
				test[0] = atoi(optarg);
				break;
			case 'b':
				test[1] = atoi(optarg);
				break;
			case 'c':
				test[2] = atoi(optarg);
				break;
			case 'd':
				test[3] = atoi(optarg);
				break;
			default:
				break;
			}
		
		
		/*cout<<"opt: "<<opt<<endl;
		cout<<"optarg: "<<optarg<<endl;
		temp = optarg;
		test[i] = atoi(optarg);
		cout<<i<<" element = "<<temp<<" "<<test[i]<<endl;
		i = i + 1;*/

	}
	int *test1 = new int[test[1]];
	test1[0] = 25;
	int b = 100000000;
	stringstream ss;
	ss.precision(1);
	ss.setf(ios::scientific);
	//ss<<std::setprecision(3)<<scientific;
	ss<<"results1\\";
	ss<<double(b);
	ss<<"-";
	ss<<test1[0];
	string test2;
	test2 = ss.str().c_str();
	cout<<" test output: "<<test2<<endl;
	cout<<test[0]<<" "<<test[1]<<" "<<test[2]<<" "<<test[3]<<endl;
    /*int a[3][2];
    int b[num1+num2] = {0,2};
    int *ptr;
    int (*p)[2];
    a[0][0] = 1;
    a[0][1] = 23;
    a[1][0] = 3;
    a[1][1] = 4;
    a[2][0] = 5;
    a[2][1] = 6;
    //int *b;
    p = a;
    ptr = b;
    
    //srand((unsigned)time(NULL));
    cout << p[1] << " " <<endl;//<< rand()%20;
    //cout << *(p[4]+1) << " " <<endl;
    testprint(p);
    //a[0][0] = a[0][0] + 1;
    //cout << a[0][0];*/
    return 0;
} 


