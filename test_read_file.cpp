#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
using namespace std;
int output[2][5];
int (*read_conf(int ind))[5];
int main ()
{
    //int (*array)[5];
    int array[2][5];
    array = read_conf(1);
    for(int i=0;i<2;i++)
    {
        for(int j=0;j<5;j++)
            cout<<array[i][j]<<"  ";
        cout<<endl;
    }
    /*fstream infile;
    int array[2][5];
    int i=0;
    int j=0;
    char cNum[256] ;
    infile.open ("mol1.txt", ifstream::in);
    if (infile.is_open())
    {
         while (infile.good())
             {
                  infile.getline(cNum, 256, ',');
                  array[j][i]= atoi(cNum) ;
                  cout<<j<<"  "<<i<<"  "<<array[j][i]<<" ";
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
    }*/
    return 0;
}

int (*read_conf(int ind))[5]
{
    fstream infile;
    int i=0;
    int j=0;
    char cNum[256] ;
    infile.open ("mol1.txt", ifstream::in);
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
