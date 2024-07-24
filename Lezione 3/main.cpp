#include "random.h"
#include <iostream>

using namespace std;

int main() {

  int throws = 10000;
  int blocks = 100;

//Definizione parametri del mio problema
  int T = 1;
  double S0 = 100;
  double K = 100;
  double r = 0.1;
  double sigma = 0.25;
  Random rand;

    //Chiamata alle funzioni
  vector<double> CallOption =
      CallOption_Dir(throws, blocks, T, S0, K, r, sigma, rand)[0];
  vector<double> CallOptionErrore =
      CallOption_Dir(throws, blocks, T, S0, K, r, sigma, rand)[3];

  vector<double> PutOption =
      PutOption_Dir(throws, blocks, T, S0, K, r, sigma, rand)[0];
  vector<double> PutOptionErrore =
      PutOption_Dir(throws, blocks, T, S0, K, r, sigma, rand)[3];

  vector<double> CallOptionDisc =
      CallOption_Disc(throws, blocks, T, S0, K, r, sigma, rand)[0];
  vector<double> CallOptionDiscErrore =
      CallOption_Disc(throws, blocks, T, S0, K, r, sigma, rand)[3];

  vector<double> PutOptionDisc =
      PutOption_Disc(throws, blocks, T, S0, K, r, sigma, rand)[0];

  vector<double> PutOptionDiscErrore =
      PutOption_Disc(throws, blocks, T, S0, K, r, sigma, rand)[3];

    //Stampo valori su file
  WriteonFile("callput.dat", CallOption_Dir(throws, blocks, T, S0, K, r, sigma, rand), PutOption_Dir(throws, blocks, T, S0, K, r, sigma, rand), CallOption_Disc(throws, blocks, T, S0, K, r, sigma, rand), PutOption_Disc(throws, blocks, T, S0, K, r, sigma, rand));

  return 0;
}