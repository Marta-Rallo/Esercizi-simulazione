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
};

//Dichiarazione delle funzioni
vector<double> MediaBlocchi_Rand(double numblocchi, double throws, Random rand);
vector<double> mediablocchi(double numblocchi, double throws,
                            vector<double> media);
vector<double> mediablocchi2(double numblocchi, double throws,
                             vector<double> dati);
vector<double> Varianza(double nblocchi, vector<double> mediablocchi,
                        vector<double> mediablocchi2);
void WriteonFile(string nomeFile, const vector<double> sumprog,
                 const vector<double> var, const vector<double> Sumprog,
                 const vector<double> Var, const vector<double> chi2);
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
