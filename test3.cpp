#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;
  
int main()
{
    int a[2][2];
    a[0][0] = 1;
    a[0][1] = 2;
    a[1][0] = 3;
    a[1][1] = 4;
    int *b;
    b = a[1];
    srand((unsigned)time(NULL));
    cout << a[0][0] << " " <<endl;//<< rand()%20;
    a[0][0] = a[0][0] + 1;
    cout << a[0][0];
    return 0;
} 
