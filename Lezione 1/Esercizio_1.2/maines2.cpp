#include "random.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main() {

  Random rnd;

  int throws = 10000;

  fstream foutput; // foutput è il nome che do alla variabile
  vector<int> N = {1, 2, 10, 100};
  vector<string> datafile = {"n1.dat", "n2.dat", "n10.dat", "n100.dat"};

  for (int i = 0; i < N.size(); i++) {
    Generate_Data(N[i], datafile[i], throws, rnd);
  }

  return 0;
}
