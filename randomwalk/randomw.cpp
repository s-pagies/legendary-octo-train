#include <fstream>
#include <iostream>
#include <vector>
#include <random>
using namespace std;

int main()
{
  std::mt19937 gen;
  gen.seed(42);

  uniform_real_distribution<double> dis(0,1);

  int anzahl=1000000;
  vector<int> ortsvektor(anzahl,0);
  int schritte=100;
  vector<float> ortsbetrag(schritte,0);
  vector<float> quadratsbetrag(schritte,0);

  for (int n=0; n<schritte; n++)
  {
    float summe = 0.0;
    float quadratsumme = 0.0;

    for (int i=0; i<anzahl; i++)
    {
      double zufall=dis((gen));
      if (zufall<0.5){
	ortsvektor[i]++;}
      else{
	ortsvektor[i]--;}
      
      summe = summe + ortsvektor[i];
      quadratsumme = quadratsumme + ortsvektor[i]*ortsvektor[i];
    }
    ortsbetrag[n] = summe/anzahl;
    quadratsbetrag[n] = quadratsumme/anzahl;
  }
  ofstream out("ortsbetrag.txt");
  for (int o=0; o<schritte; o++)
  {
    out<<ortsbetrag[0]<<endl;
  }
  ofstream arg("quadratsbetrag.txt");
  for (int o=0; o<schritte; o++)
  {
    arg<<quadratsbetrag[o]<<endl;
  }
}
