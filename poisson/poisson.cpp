#include <iostream>
#include <array>
#include <fstream>
#include <cmath>

using namespace std;

double epsilon_0  = 8.854*pow(10,-12);
double epsilon_r  = 1;
double epsilon	  = epsilon_0*epsilon_r;
double rho	  = 0;
double constant	  = rho/epsilon;
double sixth	  = 1.0/6.0;

// define 3D tensor
const int i = 80;
const int j = 80;
const int k = 80;
array<array<array<double,i+1>,j+1>,k+1> grid;

double steps  = 15000; //maximum number of iteration steps
double minimum_number_of_steps = 500; // speaks for itself
int fall      = 0; 
// fall = 0 gives a charged plate (upper boundary gridpoints)
// fall = 1 gives a single charged gridpoint a fixed charge
bool neumann  = 1;
// neumann = 0 deactivates neumann boundary condition -> same effect as a groundet boundary
// neumann = 0 activates neumann boundary condition -> boundary acts as insulator
double charge_of_particle = -1000; // charge of gridpoint for fall 1

double abort_condition = 1001;
double some_value_for_accuaracy = 0.002;
double number_of_executed_steps = 0;

int main()
{
  if(fall==0) // fall 0 = charged plate as boundry condition
  {
    // defining output
    ofstream out2d("2dprojection_plate.txt");
    ofstream out1d("1dprojection_plate.txt");
    // set uppermost gridpoits to  1000V leaving all others at 0
    for(int aa=1;aa<i+1;aa++)
    {
      for(int bb=1;bb<j+1;bb++){
        grid[1][aa][bb]=1000;}
    }

    // loop over all steps
    for(int d=1; d<steps+1; d++)
    {
      // calculate next steps
      for(int e=2;e<i+1-1;e++)
      {
	for(int f=1;f<j+1;f++)
	{
	  for(int g=1;g<k+1;g++)
	  {
	  // this is the central formula calculating the electrical potential:
	    grid[e][f][g] = sixth*((grid[e+1][f][g] + grid[e-1][f][g] + grid[e][f+1][g] + grid[e][f-1][g] + grid[e][f][g+1] + grid[e][f][g-1]) + constant);
	    if (neumann==true)
	    {
	      // applying the neumann boundary condition for each step
	      if (f==j-1){
		  grid[e][f+1][g] = grid[e][f][g];
		  continue;}
	      if (g==k-1){
		  grid[e][f][g+1] = grid[e][f][g];
		  continue;}
	      if (f==1){
		  grid[e][f-1][g] = grid[e][f][g];
		  continue;}
	      if(g==1){
		  grid[e][f][g-1] = grid[e][f][g];
		  continue;}
	    }
	  }
	}
      }
      number_of_executed_steps = d;
      // abort condition
      double wert  = abs(grid[i/2][j/2][k/2] - abort_condition);
      if (wert<some_value_for_accuaracy && d>minimum_number_of_steps){
	break;}
      abort_condition = grid[i/2][j/2][k/2];
    }// end of steps

    // write results into txt-files
    for(int aa=1;aa<i+1;aa++){
      out1d  << aa << ", " << grid[aa][j/2][k/2] << endl;}
    for(int aa=1;aa<i+1;aa++)
    {
      for(int bb=1;bb<j+1;bb++){
        out2d  << aa << ", " << bb << ", " << grid[aa][bb][k/2] << endl;}
    }

    // create format.txt for plotting the data (number of iteration steps, gridmeasurements)
    ofstream out0("format.txt");
    out0 << number_of_executed_steps <<", " << i <<", " << j <<", " << k <<endl;
    for(int e=0;e<i+1;e++)
    {
      for(int f=0;f<j+1;f++)
      {
	for(int g=0;g<k+1;g++){
	  grid[e][f][g] = 0;}
      }
    }
  }// end fall==0
  else if(fall==1)
  {
    ofstream out2d("2dprojection_point.txt");
    ofstream out1d("1dprojection_point.txt");
    // set a fixed charge for a single gridpoint
    // loop over all steps
    for(int d=1; d<steps+1; d++)
    {
      // calculate next steps
      for(int e=1;e<i+1-1;e++)
      {
	for(int f=1;f<j+1-1;f++)
	{
	  for(int g=1;g<k+1-1;g++)
	  {
	    if(e==i/2 && f==j/2 && g==k/2){
	      rho=charge_of_particle;}
	    else{rho=0;}
	      
	    grid[e][f][g] = sixth*((grid[e+1][f][g] + grid[e-1][f][g] + grid[e][f+1][g] + grid[e][f-1][g] + grid[e][f][g+1] + grid[e][f][g-1]) + constant);
	    if (neumann==true) // neumann boundary condition again now for all 6 sides of the cube
	    {
	      if (e==i-2){
		  grid[e+1][f][g] = grid[e][f][g];}
	      if (f==j+1){
		  grid[e][f+1][g] = grid[e][f][g];
		  continue;}
	      if (g==k+1){ 
		  grid[e][f][g+1] = grid[e][f][g];
		  continue;}
	      if (e==1){
		  grid[e-1][f][g] = grid[e][f][g];}
	      if (f==1){
		  grid[e][f-1][g] = grid[e][f][g];
		  continue;}
	      if(g==1){
		  grid[e][f][g-1] = grid[e][f][g];
		  continue;}
	    }
	  }
	}
      }
    }// end of steps

    // writing output into txt-files
    for(int aa=1;aa<i+1;aa++){
      out1d  << aa << ", " << grid[aa][j/2][k/2] << endl;}
    for(int aa=1;aa<i+1;aa++)
    {
      for(int bb=1;bb<j+1;bb++){
      out2d  << aa << ", " << bb << ", " << grid[aa][bb][k/2] << endl;}
    }
  }//end of fall 1
  return 0;
}
