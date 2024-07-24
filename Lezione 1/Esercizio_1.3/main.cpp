#include "EsperimentoBuffon.h"
#include "random.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char **argv) {
  // Controlla se il numero di argomenti passati è corretto
  if (argc != 3) {
    cout << "Errore! Corretto utilizzo del programma è " << argv[0]
         << " , il valore di d (distanza delle righe nel piano) e L (lunghezza "
            "dell'ago)."
         << endl;
    return -1; // Termina il programma se il numero di argomenti è errato
  }

  int points = 1000000;        // Numero totale di lanci
  int blocks = 100;            // Numero di blocchi
  int lanci = points / blocks; // Numero di lanci per blocco

  double d = atof(
      argv[1]); // Converte il primo argomento in double (distanza delle righe)
  double L = atof(
      argv[2]); // Converte il secondo argomento in double (lunghezza dell'ago)

  EsperimentoBuffon Exp(
      d, L); // Crea un oggetto EsperimentoBuffon con i parametri forniti

  vector<double> pi = Exp.Esegui(
      points, blocks); // Esegue l'esperimento e ottiene le stime di pi greco

  vector<double> pi_mediablocchi =
      mediablocchi(blocks, lanci, pi); // Calcola le medie dei blocchi
  vector<double> pi_mediablocchi2 =
      mediablocchi2(blocks, lanci, pi); // Calcola le medie quadrate dei blocchi
  vector<double> pi_err = Varianza(blocks, pi_mediablocchi,
                                   pi_mediablocchi2); // Calcola la varianza

  WriteonFile("pi.dat", pi_mediablocchi,
              pi_err); // Scrive i risultati su un file

  return 0; // Termina il programma
}
