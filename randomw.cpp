#include <fstream>
#include <iostream>
#include <vector>
#include <random>
using namespace std;

//loop: get a random number to decide to step left or right
//
int main()
{
  std::mt19937 gen;
  gen.seed(42);

  uniform_real_distribution<double> dis(0,1);

  x=0;
  int schritte=100;
  int anzahl=1;

  for (int n=0; n<schritte; ++n)
  {
    double z=dis((gen));
    if (dis((gen))<0.5)
      {x=x+1;}
    else
      {x=x-1;}
    cout<<x<<endl;
  }
}
