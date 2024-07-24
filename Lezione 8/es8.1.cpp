#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

// Define constants
const double hbar = 1.0;
const double m = 1.0;

// Trial wave function
double Psi_trial(double x, double sigma, double mu) {
    return exp(-pow(x - mu, 2) / (2 * sigma * sigma)) +
           exp(-pow(x + mu, 2) / (2 * sigma * sigma));
}

// Square of trial wave function
double Psi_squared(double x, double sigma, double mu) {
    double psi = Psi_trial(x, sigma, mu);
    return psi * psi;
}

// Energy
double Energy(double x, double sigma, double mu) {
    double sigma2 = sigma * sigma;
    double exp1 = exp(-pow(x - mu, 2) / (2 * sigma2));
    double exp2 = exp(-pow(x + mu, 2) / (2 * sigma2));
    double kinetic = (exp1/sigma2)*(((x-mu)*(x-mu)/sigma2)-1)+(exp2/sigma2)*(((x+mu)*(x+mu)/sigma2)-1);
    double potential = pow(x, 4) - 2.5 * x * x;
    return (-((hbar*hbar)/(2*m))*kinetic)/(Psi_trial(x,sigma,mu)) + potential;
}

int main() {
    // Parameters
    double sigma = 0.627525;
    double mu = -0.780447;
    int nsteps = 1000000; // Total number of Metropolis steps
    int nblocks = 10000; // Number of blocks
    Random rand;
    vector <double> dati(nblocks);

    // Metropolis parameters
    double step_size = 2.5; // Uniform step size for Metropolis algorithm

    // Initial position
    double x = 0.0;
    int accepted_moves = 0;
    // Perform Metropolis algorithm
    for (int i = 0; i < nblocks; i++) {
        double energy = 0;
        for(int j=0; j < nsteps/nblocks; j++) {
        // Propose a new position
        double x_new = x + rand.Rannyu(-step_size,step_size);

        // Calculate the acceptance probability
        double ratio = Psi_squared(x_new, sigma, mu) / Psi_squared(x, sigma, mu);
        double acceptance_prob = min(1.0, ratio);

        // Accept or reject the move
        if (rand.Rannyu() < acceptance_prob) {
            x = x_new;
            accepted_moves++;
        }

        // Calculate energy
        energy += Energy(x, sigma, mu);

        }

        dati[i] = energy/(nsteps/nblocks);

    }


    vector<vector<double>> dati_elaborati = elaboraDati(nblocks, dati);
    scriviSuFile("risultati_8.1.dat", {dati_elaborati[0], dati_elaborati[2]});
    double acceptance_rate = static_cast<double>(accepted_moves) / nsteps;
    cout << "Acceptance rate: " << acceptance_rate << endl;

    return 0;
}
