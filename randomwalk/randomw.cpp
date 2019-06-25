#include <fstream>
#include <iostream>
#include <vector>
#include <random>
using namespace std;

int main()
{
//integrieren des Zufallszahlgenerators
  std::mt19937 gen;
  gen.seed(42);

  uniform_real_distribution<double> dis(0,1);

//festlegen der Anzahlen von Schritten, Teilchen sowie der verschiedenen Vektoren
  int anzahl=1000000;
  vector<int> ortsvektor(anzahl,0);
  int schritte=100;
  vector<float> ortsbetrag(schritte,0);
  vector<float> quadratsbetrag(schritte,0);
  vector<int> dichte(schritte,0);

//Schleife über alle Schritte
  for (int n=0; n<schritte; n++)
  {
    float summe = 0.0;
    float quadratsumme = 0.0;

//Schleife über alle Teilchen
    for (int i=0; i<anzahl; i++)
    {
      double zufall=dis((gen));
      if (zufall<0.5){
	ortsvektor[i]++;}
      else{
	ortsvektor[i]--;}

//Berechnung von <x> und <x²>      
      summe = summe + ortsvektor[i];
      quadratsumme = quadratsumme + ortsvektor[i]*ortsvektor[i];
    }
    ortsbetrag[n] = summe/anzahl;
    quadratsbetrag[n] = quadratsumme/anzahl;
  }

//Dichteberechnung am Ende aller Schritte
  for (int z=0; z<schritte; z++)
  {
    for (int j=0; j<anzahl; j++)
    {
      if (ortsvektor[j]=z){
	dichte[z]++;}
    }  
  }
//Schreiben der Ergebnisse in txt-Dateien
  ofstream out("ortsbetrag.txt");
  for (int o=0; o<schritte; o++)
  {
    out<<ortsbetrag[o]<<endl;
  }
  ofstream arg("quadratsbetrag.txt");
  for (int o=0; o<schritte; o++)
  {
    arg<<quadratsbetrag[o]<<endl;
  }
  ofstream di("dichte.txt");
  for (int o=0; o<schritte; o++)
  {
    di<<dichte[o]<<endl;
  }
}
