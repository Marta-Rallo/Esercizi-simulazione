#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

int main() {
		Random rnd;

		int M = 100000; // Numero di lanci
		int nblocchi = 100; // Numero di blocchi

		// Chiamata alla funzione media e inizializzazione del vettore Media
		vector<double> Media = MediaBlocchi_Rand(nblocchi, M, rnd);

		// Chiamata alla funzione mediablocchi e inizializzazione del vettore Mediablocchi
		vector<double> Mediablocchi = mediablocchi(nblocchi, M, Media);

		// Chiamata alla funzione mediablocchi2 e inizializzazione del vettore Mediablocchi2
		vector<double> Mediablocchi2 = mediablocchi2(nblocchi, M, Media);

		// Chiamata alla funzione errore e inizializzazione del vettore Err
		vector<double> Err = Varianza(nblocchi, Mediablocchi, Mediablocchi2);

		double L = M / nblocchi; // Numero di lanci in ogni blocco, assicurati che M sia un multiplo di nblocchi
		double num = 0; //Variabile d'appoggio

		// Vettore per le medie dei blocchi 
		vector<double> ave(nblocchi);

		// Calcola le medie dei blocchi
		for (int k = 0; k < nblocchi; k++) {
				num = 0;

				for (int i = 0; i < L; i++) {
						// Calcola il quadrato della differenza tra il numero generato casualmente e 0.5 e sommalo
						num += pow(rnd.Rannyu() - 0.5, 2);
				}   
				// Calcola la media del blocco k
				ave[k] = num / L;
		}

		// Calcola le medie cumulative dei blocchi e i quadrati delle medie cumulative dei blocchi
		vector<double> Sumprog = mediablocchi(nblocchi, M, ave);
		vector<double> Sumprog2 = mediablocchi2(nblocchi, M, ave);

		// Calcola l'errore utilizzando le medie cumulative dei blocchi e i quadrati delle medie cumulative dei blocchi
		vector<double> Var = Varianza(nblocchi, Sumprog, Sumprog2);

		//Test chi2
		vector<double> vec_chi2 = chi2(100, 10000, 100, rnd);

		// Scrivi i risultati su file
		WriteonFile("ris.dat", Mediablocchi, Err, Sumprog, Var, vec_chi2);



		return 0;
}