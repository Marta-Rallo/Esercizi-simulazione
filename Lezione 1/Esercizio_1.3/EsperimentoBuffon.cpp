#include "EsperimentoBuffon.h"
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

// Controlla i parametri dell'esperimento per assicurarsi che siano validi
void EsperimentoBuffon::ControllaParametri(double dist, double lunghago) {
    if (dist < lunghago || dist >= pow(lunghago, 10)) {
        cout << "Errore! La distanza fra le righe del piano deve essere maggiore "
                "della lunghezza del'ago, ma non oltre la lunghezza dell'ago alla "
                "decima!"
             << endl;
        exit(0);
    }
}

// costruttore con dati che passerò da main
EsperimentoBuffon::EsperimentoBuffon(double dist, double lunghago) : rdn() {
    ControllaParametri(dist, lunghago); // Controlla i parametri passati

    d = dist;    // Imposta la distanza tra le righe del piano
    L = lunghago; // Imposta la lunghezza dell'ago
}

// Funzione che esegue l'esperimento di Buffon
vector<double> EsperimentoBuffon::Esegui(int throws, int blocks) {
    double xi, xf;
    double lanci = float(throws / blocks); // lanci in ogni blocco
    double L = GetLAgo(); // Ottiene la lunghezza dell'ago
    double d = GetDPlane(); // Ottiene la distanza tra le righe del piano
    double estimate_pi; // Variabile per stimare il valore di pi greco
    vector<double> pi(blocks); // Vettore per memorizzare le stime di pi greco per ogni blocco

    // per simmetrie del sistema mi riduco a una striscia verticale
    for (int i = 0; i < blocks; i++) {
        int hit = 0; // Contatore per il numero di volte in cui l'ago colpisce una riga

        // Esegue i lanci in ogni blocco
        for (int j = 0; j < lanci; j++) {
            xi = rdn.Rannyu(0, d); // Genera una posizione iniziale casuale per l'ago
            double cos = rdn.Cos(); // Genera un angolo casuale (coseno dell'angolo)
            xf = xi + L * cos; // Calcola la posizione finale dell'ago

            // Verifica se l'ago colpisce una riga
            if (xf > d || xf < 0) {
                hit++;
            }
        }

        // Stima il valore di pi greco per il blocco corrente
        pi[i] = estimate_pi = (2 * L * lanci) / (hit * d);
    }

    return pi; // Ritorna il vettore con le stime di pi greco per ogni blocco
}

