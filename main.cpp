#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

double fahrenheit(double celsius){
  return celsius*9./5.+32;
}

//double celsius(double fahrenheit){
//  return fahrenheit....

int main(){

  //liest input.txt in inputvector ein
  vector<double> inputvector;
  std::ifstream ga("input.txt");
  double var1;
  for(int i=0; i<3; i++){
    ga>>var1;
    inputvector.push_back(var1);
  }
  
  //wendet Celsius-Fahrenheit Formel auf inputvector an und schreibt Ergebnis in outputvector
  vector<double> outputvector;
  double var2;
  for(int i=0; i<inputvector.size(); i++){
    double var2 = fahrenheit(inputvector[i]);
    outputvector.push_back(var2);
  }

  //writes outputvector into output.txt on harddrive
  ofstream out("output.txt");
  for(int i=0; i<outputvector.size(); i++){
    out << outputvector[i] << endl;
  }

  return 0; 
}
