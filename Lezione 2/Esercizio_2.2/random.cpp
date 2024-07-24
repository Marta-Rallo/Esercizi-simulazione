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

double Random ::Rannyu(double min, double max) {
  // This function generates a random number in the range [min, max)
  return min + (max - min) * Rannyu();
}

int Random ::Rannyu_Int(int min, int max) {

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

// Esegue Random Walk nel discreto
vector<vector<double>> MediaPosition_Rand(int throws, int blocks, int steps,
                                          Random rand, double a) {

  vector<vector<double>> MeanBlocks;
  vector<double> mediaposition(steps);
  vector<double> mediaposition2(steps);
  vector<double> mediapositionsquare(steps);
  vector<double> mediapositionerr(steps);

  double sumxf, sumyf, sumzf;
  double dist;
  double app;
  int L = throws / blocks;

  for (int passo = 1; passo <= steps; passo++) {

    for (int i = 0; i < blocks; i++) {
      app = 0;

      for (int j = 0; j < L; j++) {
        sumxf = 0;
        sumyf = 0;
        sumzf = 0;

        for (int k = 0; k < passo; k++) {

          // scelgo la direzione in cui muovermi fra 1, 2, 3 ovvero x, y, z
          // rispettivamente: per isotropia dello spazio scelgo la direzione
          // secondo una distribuzione uniforme

          int direction = rand.Rannyu_Int(1, 4);

          if (direction == 1) {
            int step =
                rand.Rannyu_Int(1, 3); // associo a 1: vado avanti di +1,
                                       // associo a 2: vado indietro di -1
            if (step == 1) {
              sumxf += a;
            } else {
              sumxf += -a;
            }

          } else if (direction == 2) {
            int step = rand.Rannyu_Int(1, 3);
            if (step == 1) {
              sumyf += a;
            } else {
              sumyf += -a;
            }

          } else {
            int step = rand.Rannyu_Int(1, 3);
            if (step == 1) {
              sumzf += a;
            } else {
              sumzf += -a;
            }
          }
        }

        dist = sqrt(pow(sumxf, 2) + pow(sumyf, 2) + pow(sumzf, 2));

        app += dist;
      }
      // Calcola media nel blocco attuale
      app /= L;
      mediaposition[passo - 1] += app;
      mediaposition2[passo - 1] += app * app;
    }

    mediaposition[passo - 1] = mediaposition[passo - 1] / steps;
    mediaposition2[passo - 1] = mediaposition2[passo - 1] / steps;
    mediapositionsquare[passo - 1] =
        mediaposition[passo - 1] * mediaposition[passo - 1];

    if (passo == 1) {
      mediapositionerr[passo - 1] = 0;
    } else {
      mediapositionerr[passo - 1] =
          sqrt((mediaposition2[passo - 1] - mediapositionsquare[passo - 1]) /
               (steps - 1));
    }
  }

  // Aggiunta di un vettore al vettore totale MeanBlocks
  MeanBlocks.push_back(mediaposition);
  MeanBlocks.push_back(mediaposition2);
  MeanBlocks.push_back(mediapositionsquare);
  MeanBlocks.push_back(mediapositionerr);

  return MeanBlocks;
}

// Eseguo Random Walk nel continuo
vector<vector<double>> MediaPosition_Cont(int throws, int blocks, int steps,
                                          Random rand, double a) {

  vector<vector<double>> MeanBlocks;
  vector<double> mediaposition(steps);
  vector<double> mediaposition2(steps);
  vector<double> mediapositionsquare(steps);
  vector<double> mediapositionerr(steps);

  double sumxf, sumyf, sumzf;
  int L = throws / blocks;
  double app;
  double dist;

  for (int passo = 1; passo <= steps; passo++) {

    for (int i = 0; i < blocks; i++) {
      app = 0;

      for (int j = 0; j < L; j++) {
        sumxf = 0;
        sumyf = 0;
        sumzf = 0;

        for (int k = 0; k < passo; k++) {
          double phi, sintheta, costheta, theta;

          // genero angolo phi uniformemenete distribuito in [0,2pi]
          phi = rand.Rannyu(0, 2 * M_PI);
          // genero sentheta uniformemenete distribuito in [0,pi]

          costheta = rand.Rannyu(-1, 1);
          sintheta = sqrt(1 - pow(costheta, 2));

          sumxf += a * cos(phi) * sintheta;
          sumyf += a * sin(phi) * sintheta;
          sumzf += a * costheta;
        }

        dist = sqrt(pow(sumxf, 2) + pow(sumyf, 2) + pow(sumzf, 2));

        app += dist;
      }

      // Calcola media nel blocco attuale
      app /= L;

      mediaposition[passo - 1] += app;
      mediaposition2[passo - 1] += app * app;
    }

    mediaposition[passo - 1] = mediaposition[passo - 1] / blocks;
    mediaposition2[passo - 1] = mediaposition2[passo - 1] / blocks;
    mediapositionsquare[passo - 1] =
        mediaposition[passo - 1] * mediaposition[passo - 1];

    if (passo == 1) {
      mediapositionerr[passo - 1] = 0;
    } else {
      mediapositionerr[passo - 1] =
          sqrt((mediaposition2[passo - 1] - mediapositionsquare[passo - 1]) /
               (steps - 1));
    }
  }

  // Aggiunta di un vettore al vettore totale MeanBlocks
  MeanBlocks.push_back(mediaposition);
  MeanBlocks.push_back(mediaposition2);
  MeanBlocks.push_back(mediapositionsquare);
  MeanBlocks.push_back(mediapositionerr);

  return MeanBlocks;
}

void WriteonFile(string nomeFile, const vector<double> sumprog,
                 const vector<double> Var, const vector<double> mean2,
                 const vector<double> err2) {

  fstream foutput;                  // foutput è il nome che do alla variabile
  foutput.open(nomeFile, ios::out); // risultato.dat è il nome del file
  if (foutput.good()) // controllo che funzioni e gli faccio stampare,
                      // altrimenti gli faccio dare errore
  {
    foutput << "MBlocks_Disc"
            << " "
            << "EBlocks_Disc"
            << " "
            << "MBlocks_Cont"
            << " "
            << "EBlocks_Cont" << endl;

    for (double i = 0; i < sumprog.size(); i++) {
      foutput << sumprog[i] << " " << Var[i] << " " << mean2[i] << " "
              << err2[i] << endl;
    }
  } else
    cout << "Errore: file is not good!" << endl;

  foutput.close(); // chiudo il file
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
