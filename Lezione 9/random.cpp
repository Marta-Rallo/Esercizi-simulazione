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

double Random ::CauchyLor(double gamma, double mu) {

	double x = Rannyu();
	double y;
	if (x == 0 || x == 1) {
		y = CauchyLor(gamma, mu);
	} else {
		y = (mu + gamma * tan(M_PI * (x - (1. / 2.))));
	}
	return y;
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

vector<double> MediaBlocchi_Rand(double numblocchi, double throws,
																 Random rand) {
	double L = throws / numblocchi; // Number of throws in each block, please use
																	// for M a multiple of N
	double sum = 0;                 // somme parziali
	vector<double> vec_media(numblocchi);

	for (int k = 0; k < numblocchi; k++) {
		sum = 0;

		for (int i = 0; i < L; i++) {
			sum += rand.Rannyu();
		}
		vec_media[k] = sum / L;
	}

	return vec_media;
}

vector<double> MediaBlocchi_minmax(double numblocchi, double throws,
																	 Random rand, double min, double max) {
	double L = throws / numblocchi; // Number of throws in each block, please use
																	// for M a multiple of N
	double sum = 0;                 // somme parziali
	vector<double> vec_media(numblocchi);

	for (int k = 0; k < numblocchi; k++) {
		sum = 0;

		for (int i = 0; i < L; i++) {
			sum += rand.Rannyu(min, max);
		}
		vec_media[k] = sum / L;
	}

	return vec_media;
}

vector<double> MediaBlocchi_Exp(double numblocchi, double throws, Random rand,
																double mean) {
	double L = throws / numblocchi; // Number of throws in each block, please use
																	// for M a multiple of N
	double sum = 0;                 // somme parziali
	vector<double> vec_media(numblocchi);

	for (int k = 0; k < numblocchi; k++) {
		sum = 0;

		for (int i = 0; i < L; i++) {
			sum += rand.Exp(mean);
		}
		vec_media[k] = sum / L;
	}

	return vec_media;
}

vector<double> mediaCL(double numblocchi, double throws, Random rand,
											 double gamma, double mu) {
	double L = throws / numblocchi; // Number of throws in each block, please use
																	// for M a multiple of N
	double sum = 0;                 // somme parziali
	vector<double> vec_media(numblocchi);

	for (int k = 0; k < numblocchi; k++) {
		sum = 0;

		for (int i = 0; i < L; i++) {
			sum += rand.CauchyLor(gamma, mu);
		}
		vec_media[k] = sum / L;
	}

	return vec_media;
}

vector<double> mediablocchi(double numblocchi, double throws,
														const vector<double> media) {

	vector<double> sumprog(numblocchi);

	for (int i = 0; i < numblocchi; i++) {
		for (int j = 0; j <= i; j++) {

			sumprog[i] += media[j];
		}
		sumprog[i] /= (i + 1);
	}

	return sumprog;
}

vector<double> mediablocchi2(double numblocchi, double throws,
														 const vector<double> dati) {

	vector<double> vec_media2(numblocchi);
	vector<double> sumprog2(numblocchi);

	for (int k = 0; k < numblocchi; k++) {
		vec_media2[k] = dati[k] * dati[k];
	}

	for (int i = 0; i < numblocchi; i++) {
		for (int j = 0; j <= i; j++) {

			sumprog2[i] += vec_media2[j];
		}
		sumprog2[i] /= (i + 1);
	}

	return sumprog2;
}

vector<double> Varianza(double nblocchi, vector<double> mediablocchi,
												vector<double> mediablocchi2) {

	vector<double> var(nblocchi);
	for (int i = 0; i < nblocchi; i++) {

		if (i == 0) {
			var[i] = 0;
		} else {
			var[i] = sqrt((mediablocchi2[i] - pow(mediablocchi[i], 2)) / (i));
		}
	}

	return var;
}

void Generate_Data(int N, string datafile, int throws, Random rnd) {

	ofstream outfile;
	double sum_rannyu;
	double sum_exp;
	double sum_cl;

	outfile.open(datafile);

	if (outfile.good()) {

		outfile << "Unif"
						<< " "
						<< "Exp"
						<< " "
						<< "CLorentz" << endl;

		for (int i = 0; i < throws; i++) {

			sum_rannyu = 0;
			sum_exp = 0;
			sum_cl = 0;

			for (int j = 1; j <= N; j++) {
				sum_rannyu += rnd.Rannyu();
				sum_exp += rnd.Exp(1);
				sum_cl += rnd.CauchyLor(1, 0);
			}

			outfile << sum_rannyu / N << " " << sum_exp / N << " " << sum_cl / N
							<< endl;
		}
	} else {
		cout << "Errore: file is not good!" << endl;
	}

	outfile.close();
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
