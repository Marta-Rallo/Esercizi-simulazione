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
#include <algorithm>
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

double Random ::Gauss(double mean, double sigma) {
  // This function generates a random number from a Gaussian distribution with
  // given mean and sigma
  double s = Rannyu();
  double t = Rannyu();
  double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
  return mean + x * sigma;
}

// Calcola l'opzione Call campionando direttamente il prezzo finale dell'asset
// S(t)
vector<vector<double>> CallOption_Dir(int throws, int blocks, int T, double S0,
                                      double prezzopattuito, double r,
                                      double sigma, Random rand) {

  vector<vector<double>> MeanBlocks;
  vector<double> Calloption(blocks);
  vector<double> Calloptionsquare(blocks);
  vector<double> Calloption2(blocks);
  vector<double> Errore(blocks);

  int L = throws / blocks;
  double Call_Dir;
  double sum;
  double sum2;
  double Call_Dirsquare;

  for (int i = 0; i < blocks; i++) {
    sum = 0;
    sum2 = 0;
    for (int j = 0; j < L; j++) {

      double Wt = rand.Gauss(0, T);
      double St = S0 * exp((r - (0.5 * sigma * sigma)) * T + (sigma * Wt));

      double profit_at_expiry = max(0.0, St - prezzopattuito);

      Call_Dir = exp(-r * T) * profit_at_expiry;
      sum = sum + Call_Dir;

      Call_Dirsquare = Call_Dir * Call_Dir;
      sum2 = sum2 + Call_Dirsquare;
    }

    Calloption[i] = sum / L;
    Calloption2[i] = Calloption[i] * Calloption[i];
    Calloptionsquare[i] = sum2 / L;

    if (i == 0) {
      Errore[i] = 0;
    } else {
      Errore[i] = sqrt((Calloptionsquare[i] - Calloption2[i]) / (blocks - 1));
    }
  }

  // Aggiunta di un vettore al vettore totale MeanBlocks
  MeanBlocks.push_back(Calloption);
  MeanBlocks.push_back(Calloptionsquare);
  MeanBlocks.push_back(Calloption2);
  MeanBlocks.push_back(Errore);

  return MeanBlocks;
}

// Calcola l'opzione Put campionando direttamente il prezzo finale dell'asset
// S(t)
vector<vector<double>> PutOption_Dir(int throws, int blocks, int T, double S0,
                                     double prezzopattuito, double r,
                                     double sigma, Random rand) {

  vector<vector<double>> MeanBlocks;
  vector<double> Putoption(blocks);
  vector<double> Putoptionsquare(blocks);
  vector<double> Putoption2(blocks);
  vector<double> Errore(blocks);

  int L = throws / blocks;
  double Put_Dir;
  double Put_Dirsquare;
  double sum;
  double sum2;

  for (int i = 0; i < blocks; i++) {
    sum = 0.;
    sum2 = 0.;
    for (int j = 0; j < L; j++) {

      double Wt = rand.Gauss(0., T);
      double St = S0 * exp((r - (0.5 * sigma * sigma)) * T + (sigma * Wt));

      double profit_at_expiry = max(0.0, prezzopattuito - St);

      Put_Dir = exp(-r * T) * profit_at_expiry;
      sum = sum + Put_Dir;

      Put_Dirsquare = Put_Dir * Put_Dir;
      sum2 = sum2 + Put_Dirsquare;
    }

    Putoption[i] = sum / L;
    Putoption2[i] = Putoption[i] * Putoption[i];
    Putoptionsquare[i] = sum2 / L;

    if (i == 0) {
      Errore[i] = 0;
    } else {
      Errore[i] = sqrt((Putoptionsquare[i] - Putoption2[i]) / (blocks - 1));
    }
  }

  // Aggiunta di un vettore al vettore totale MeanBlocks
  MeanBlocks.push_back(Putoption);
  MeanBlocks.push_back(Putoptionsquare);
  MeanBlocks.push_back(Putoption2);
  MeanBlocks.push_back(Errore);

  return MeanBlocks;
}

// Calcola l'opzione Call campionando il percorso discretizzato del prezzo
// dell'asset
vector<vector<double>> CallOption_Disc(int throws, int blocks, int T, double S0,
                                       double prezzopattuito, double r,
                                       double sigma, Random rand) {

  vector<vector<double>> MeanBlocks;
  vector<double> CallDisc(blocks);
  vector<double> CallDisc2(blocks);
  vector<double> CallDiscquare(blocks);
  vector<double> Errore(blocks);

  int L = throws / blocks;
  double sum;
  double St_npiuuno;
  double St_n;
  double Call_Disc;
  double Call_Discquare;
  double sum2;

  for (int i = 0; i < blocks; i++) {
    sum = 0;
    sum2 = 0;
    for (int j = 0; j < L; j++) {

      // divido intervallo (0,T) in 100 sottointervallini
      for (double t = 0; t < T; t += 0.01) {
        double Wt = rand.Gauss(0, T);

        if (t == 0) {
          St_n = S0;
        } else {
          St_n = St_npiuuno;
        }

        St_npiuuno = St_n * exp((r - (0.5 * sigma * sigma)) * (0.01) +
                                (sigma * Wt) * sqrt(0.01));

        double profit_at_expiry = max(0.0, St_npiuuno - prezzopattuito);

        Call_Disc = exp(-r * (t + 0.01)) * profit_at_expiry;
        Call_Discquare = Call_Disc * Call_Disc;
      }

      sum = sum + Call_Disc;
      sum2 = sum2 + Call_Discquare;
    }

    CallDisc[i] = sum / L;
    CallDisc2[i] = sum2 / L;
    CallDiscquare[i] = CallDisc[i] * CallDisc[i];

    if (i == 0) {
      Errore[i] = 0;
    } else {
      Errore[i] = sqrt((CallDisc2[i] - CallDiscquare[i]) / (blocks - 1));
    }
  }

  // Aggiunta di un vettore al vettore totale MeanBlocks
  MeanBlocks.push_back(CallDisc);
  MeanBlocks.push_back(CallDisc2);
  MeanBlocks.push_back(CallDiscquare);
  MeanBlocks.push_back(Errore);

  return MeanBlocks;
}

// Calcola l'opzione Call campionando il percorso discretizzato del prezzo
// dell'asset
vector<vector<double>> PutOption_Disc(int throws, int blocks, int T, double S0,
                                      double prezzopattuito, double r,
                                      double sigma, Random rand) {

  vector<vector<double>> MeanBlocks;
  vector<double> PutDisc(blocks);
  vector<double> PutDisc2(blocks);
  vector<double> PutDiscquare(blocks);
  vector<double> Errore(blocks);

  int L = throws / blocks;
  double sum;
  double sum2;
  double St_npiuuno;
  double St_n;
  double Put_Disc;
  double Put_Discquare;
  double profit_at_expiry;

  for (int i = 0; i < blocks; i++) {
    sum = 0;
    sum2 = 0;
    for (int j = 0; j < L; j++) {

      // divido intervallo (0,T) in 100 sottointervallini
      for (double t = 0; t < T; t += 0.01) {
        double Wt = rand.Gauss(0, T);

        if (t == 0) {
          St_n = S0;
        } else {
          St_n = St_npiuuno;
        }

        St_npiuuno = St_n * exp((r - (0.5 * sigma * sigma)) * (0.01) +
                                (sigma * Wt) * sqrt(0.01));

        profit_at_expiry = max(0., prezzopattuito - St_npiuuno);

        Put_Disc = (exp(-r * (t + 0.01))) * profit_at_expiry;
        Put_Discquare = Put_Disc * Put_Disc;
      }

      sum = sum + Put_Disc;
      sum2 = sum2 + Put_Discquare;
    }

    PutDisc[i] = sum / L;
    PutDisc2[i] = sum2 / L;
    PutDiscquare[i] = PutDisc[i] * PutDisc[i];

    if (i == 0) {
      Errore[i] = 0;
    } else {
      Errore[i] = sqrt((PutDisc2[i] - PutDiscquare[i]) / (blocks - 1));
    }
  }

  // Aggiunta di un vettore al vettore totale MeanBlocks
  MeanBlocks.push_back(PutDisc);
  MeanBlocks.push_back(PutDisc2);
  MeanBlocks.push_back(PutDiscquare);
  MeanBlocks.push_back(Errore);

  return MeanBlocks;
}

void WriteonFile(string nomeFile, vector<vector<double>> uno,
                 vector<vector<double>> due, vector<vector<double>> tre,
                 vector<vector<double>> quattro) {

  ofstream foutput;       // foutput è il nome che do alla variabile
  foutput.open(nomeFile); // risultato.dat è il nome del file
  if (foutput.good())     // controllo che funzioni e gli faccio stampare,
                          // altrimenti gli faccio dare errore
  {
    foutput << "Call1"
            << " "
            << "ECall1"
            << " "
            << "Put1"
            << " "
            << "EPut1"
            << " "
            << "Call2"
            << " "
            << "ECall2"
            << " "
            << "Put2"
            << " "
            << "EPut2" << endl;
    for (double i = 0; i < uno[0].size(); i++) {
      foutput << uno[0][i] << " " << uno[3][i] << " " << due[0][i] << " "
              << due[3][i] << " " << tre[0][i] << " " << tre[3][i] << " "
              << quattro[0][i] << " " << quattro[3][i] << endl;
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
