/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

// Random :: Random(){}
//  Default constructor, does not perform any action

Random ::~Random() {}
// Default destructor, does not perform any action

void Random ::SaveSeed() {
  // This function saves the current state of the random number generator to a
  // file "seed.out"
  ofstream WriteSeed;
  WriteSeed.open("seed.out");
  if (WriteSeed.is_open()) {
    WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4
              << endl;
    ;
  } else
    cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

// mio costruttore random
Random ::Random() {
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  SaveSeed();
}

double Random ::Gauss(double mean, double sigma) {
  // This function generates a random number from a Gaussian distribution with
  // given mean and sigma
  double s = Rannyu();
  double t = Rannyu();
  double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
  return mean + x * sigma;
}

double Random ::Rannyu(double min, double max) {
  // This function generates a random number in the range [min, max)
  return min + (max - min) * Rannyu();
}

double Random ::Rannyu(void) {
  // This function generates a random number in the range [0,1)
  const double twom12 = 0.000244140625;
  int i1, i2, i3, i4;
  double r;

  i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
  i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
  i3 = l3 * m4 + l4 * m3 + n3;
  i4 = l4 * m4 + n4;
  l4 = i4 % 4096;
  i3 = i3 + i4 / 4096;
  l3 = i3 % 4096;
  i2 = i2 + i3 / 4096;
  l2 = i2 % 4096;
  l1 = (i1 + i2 / 4096) % 4096;
  r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

  return r;
}

double Random ::Exp(double mean) {
  double y = Rannyu();
  return -(1 / mean) * log(1. - y);
}

double Random ::Cos() {
  double x = Rannyu(0, pow(10, 5));
  return cos(x);
}

double Random ::Sin() {
  double x = Rannyu(0, M_PI);
  return sin(x);
}

double Random ::CauchyLor(double gamma, double mu) {
  double x = Rannyu();

  return mu + gamma * tan(M_PI * (x - 0.5));
}

void Random ::SetRandom(int *s, int p1, int p2) {
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

vector<double> mediarannyu(double numblocchi, double throws, Random rand) {
  double L = throws / numblocchi;  // Number of throws in each block, please use
                                   // for M a multiple of N
  double sum = 0;                  // somme parziali
  vector<double> vec_media(numblocchi);
  vector<double> sumprog(numblocchi);

  for (int k = 0; k < numblocchi; k++) {
    sum = 0;

    for (int i = 0; i < L; i++) {
      sum += rand.Rannyu();
    }
    vec_media[k] = sum / L;
  }

  return vec_media;
}

vector<double> integraleave(double xmin, double xmax, int throws, int blocchi,Random &rnd) {
  double sum = 0;
  double L = throws / blocchi;
  double x = 0;
  vector<double> media_int(blocchi);

  for (int k = 0; k < blocchi; k++) {
    sum = 0;

    for (int i = 0; i < L; i++) {
      x = rnd.Rannyu(xmin, xmax);
      sum += (M_PI / 2) * cos((M_PI * x) / 2);
    }
    media_int[k] = (xmax - xmin) * (sum / L);
  }

  return media_int;
};

vector<double> integralehom(double xmin, double xmax, double fmax, Random &rnd,int throws, int blocchi) {
  double x = 0;

  double y = 0;

  double gx = 0;

  double L = throws / blocchi;

  vector<double> media_int(blocchi);

  for (int i = 0; i < blocchi; i++) {
    gx = 0;

    for (int j = 0; j < L; j++) {
      x = rnd.Rannyu(xmin, xmax);

      y = rnd.Rannyu(0, fmax);

      if (y < ((3 / M_PI) * (1 - (x * x)))) {
        gx += (M_PI / 3) * ((cos((M_PI * x) / 2)) / (1 - (x * x)));

      } else {
        j--;
      }
    }

    media_int[i] = gx / L;
  }

  return media_int;
};

void scriviSuFile(string nomeFile, const vector<vector<double>>& vettori) {
    fstream foutput(nomeFile, ios::out);

    if (foutput.good()) {
        // Determina la dimensione dei vettori
        int size = vettori[0].size();

        // Itera attraverso gli elementi dei vettori
        for (int i = 0; i < size; i++) {
            // Scrivi i valori di tutti i vettori nella stessa riga del file
            foutput << i+1 << " ";
            for (const auto& vettore : vettori) {
                foutput << vettore[i] << " ";
            }
            foutput << endl;
        }
    } else {
        cout << "Errore: impossibile aprire il file per la scrittura." << endl;
    }

    foutput.close(); // Chiudi il file
}

vector<vector<double>> elaboraDati(double numblocchi, const vector<double>& media) {

    vector<vector<double>> risultato(3, vector<double>(numblocchi));
    vector<double> sumprog(numblocchi);
    vector<double> vec_media2(numblocchi);
    
    
    for(int k = 0; k < numblocchi; k++) {
        vec_media2[k] = media[k] * media[k];
    }

    for(int i = 0; i < numblocchi; i++) {
        for (int j = 0; j <= i; j++) {
            sumprog[i] += media[j];
            risultato[1][i] += vec_media2[j];
        }
        risultato[0][i] = sumprog[i] / (i + 1); // Media a blocchi
        risultato[1][i] /= (i + 1); // Media a blocchi quadra
        if (i > 0)
            risultato[2][i] = sqrt((risultato[1][i] - risultato[0][i]*risultato[0][i]) / i); // Errore
        else
            risultato[2][i] = 0; // Errore = 0 se i = 0
    }

    return risultato;
}




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
