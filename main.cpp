#include <fstream>
#include <iostream>
using namespace std;

double fahrenheit(double celsius)
{
  return celsius*9./5.+32;
}

int main()
{
  double var1;
  ifstream in("input.txt");
  in>>var1;
  double Wert = fahrenheit(var1);
  cout << Wert << endl;
  return 0;
}
