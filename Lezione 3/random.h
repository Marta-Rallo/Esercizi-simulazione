/****************************************************************
*****************************************************************
                _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
         _/_/  _/ _/       _/       Physics Department
        _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

// This class contains functions for generating random numbers using the RANNYU
// algorithm
class Random {

private:
  int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;

protected:
public:
  // Default constructor
  // Random();
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int *, int, int);
  // Method to save the seed to a file
  void SaveSeed();
  // Metodo per generare un numero casuale tra 0 e 1, costruttore personale
  Random();
  // Funzione che genera un numero casuale tra 0 e 1
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);
};

// Dichiarazione delle funzioni; nelle CallOption e nelle PutOption sono già
// inclusi i calcoli delle Medie a Blocchi usando il metodo del Data Blocking
vector<vector<double>> CallOption_Dir(int throws, int blocks, int T, double S0,
                                      double prezzopattuito, double r,
                                      double sigma, Random rand);

vector<vector<double>> PutOption_Dir(int throws, int blocks, int T, double S0,
                                     double prezzopattuito, double r,
                                     double sigma, Random rand);

vector<vector<double>> CallOption_Disc(int throws, int blocks, int T, double S0,
                                       double prezzopattuito, double r,
                                       double sigma, Random rand);

vector<vector<double>> PutOption_Disc(int throws, int blocks, int T, double S0,
                                      double prezzopattuito, double r,
                                      double sigma, Random rand);

// Scrivi su file i dati
void WriteonFile(string nomeFile, vector<vector<double>> uno,
                 vector<vector<double>> due, vector<vector<double>> tre,
                 vector<vector<double>> quattro);
vector<double> chi2(double interv, double n, double num_iterations, Random rnd);

#endif // __Random__

/****************************************************************
*****************************************************************
                _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
         _/_/  _/ _/       _/       Physics Department
        _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
