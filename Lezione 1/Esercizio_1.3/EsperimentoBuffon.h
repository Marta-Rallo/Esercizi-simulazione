#ifndef _EsperimentoBuffon_
#define _EsperimentoBuffon_

#include "random.h"
#include <iostream>
#include <vector>
using namespace std;

class EsperimentoBuffon {

public:
  // costruttore e distruttore
  EsperimentoBuffon(double dist, double lunghago);
  ~EsperimentoBuffon() { ; };

//Dichiarazione delle funzioni
  void GetPlane(double lungh, double alt);
  double GetDPlane() { return d; }
  double GetLAgo() { return L; }
  vector <double> Esegui(int throws, int blocks);
  void ControllaParametri(double dist, double lunghago);

private:
  // generatori di numeri casuali
  Random rdn;
  double d; // distanza tra le righe del piano

  // parametri ago
  double L; // lunghezza ago, ricordiamoci dovrà essere d>L
};

#endif