#include <fstream>
#include <iostream>
#include <vector>
#include <random>
using namespace std;

int main()
{
//integration of random number generator
  mt19937 gen;
  gen.seed(45);

  uniform_real_distribution<double> dis(0,1);

//predefining variables and vecors
  int quantity=1000;
  int steps=100;
  vector<int> x_position(quantity,0);
  vector<int> y_position(quantity,0);

//loop over all steps
  for (int t=0; t<steps; t++)
  {
//loop over all particles moving particles for each step
    for (int n=0; n<quantity; n++)
    {
      double zufall=dis(gen);
      if (zufall < 0.25){
	x_position[n]++;}
      if (zufall >= 0.25 && zufall < 0.5){
	x_position[n]--;}
      if (zufall >= 0.5 && zufall < 0.75){
	y_position[n]++;}
      if (zufall >= 0.75){
	y_position[n]--;}
    }
  }
//write results into txt-files
  ofstream out("ort.txt");
  for (int i=0; i<quantity; i++)
  {
    out<<x_position[i]<<"\t"<<y_position[i]<<endl;
  }
}
