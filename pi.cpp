#include <random>
#include <iostream>
#include <chrono>
using namespace std;

//get two random numbers
int main()
{
  mt19937 gen;
//  gen.seed(chrono::high_resolution_clock::now()
//    .time_since_epoch().count());
  gen.seed(42);

  uniform_real_distribution<double> dis(0, 1);

  int j=0;
  int k=10000;

  for (int n = 0; n<k; ++n)
  {
    double x = dis((gen));
    double y = dis((gen));
    if (x*x+y*y<=1)
      { ++j;}
  }

  double pi = 4*j/k;
  cout<<pi<<endl;
}
