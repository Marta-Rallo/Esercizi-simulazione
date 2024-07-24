#include "random.h"
#include <iostream>

using namespace std;

int main() {
  Random rand;
  int throws = 10000; // numero di lanci
  int blocks = 100;   // numero di blocchi
  int steps = 100;    // numero di steps
  double a = 1;       // passo del random walk

  // Random walk discreto
  vector<double> position_discrete =
      MediaPosition_Rand(throws, blocks, steps, rand, a)[0];
  vector<double> Err = MediaPosition_Rand(throws, blocks, steps, rand, a)[3];

  // Random walk continuo
  vector<double> position_continuum =
      MediaPosition_Cont(throws, blocks, steps, rand, a)[0];
  vector<double> Errcont =
      MediaPosition_Cont(throws, blocks, steps, rand, a)[3];

  // Scrivo dati su file
  WriteonFile("pos.dat", position_discrete, Err, position_continuum, Errcont);

  return 0;
}