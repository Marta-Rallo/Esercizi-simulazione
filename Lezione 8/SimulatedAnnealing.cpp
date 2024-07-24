#include "SimulatedAnnealing.h"

SimulatedAnnealing::SimulatedAnnealing(int SA_steps, double initial_temperature, double cooling_rate, double temp_min, Random& rand)
    : SA_steps(SA_steps), initial_temperature(initial_temperature), cooling_rate(cooling_rate), temp_min(temp_min), rand(rand) {}

std::pair<double, double> SimulatedAnnealing::perform(double sigma, double mu) {
    int nsteps = 10000;
    int nblocks = 1000;
    double step_size = 2.5;
    std::vector<double> dati(nblocks);
    double x = 0.0;
    int accepted_moves = 0;

    WaveFunction wf(sigma, mu);

    for (int i = 0; i < nblocks; i++) {
        double energy = 0;
        for (int j = 0; j < nsteps / nblocks; j++) {
            double x_new = x + rand.Rannyu(-step_size, step_size);
            double ratio = wf.PsiSquared(x_new) / wf.PsiSquared(x);
            double acceptance_prob = std::min(1.0, ratio);

            if (rand.Rannyu() < acceptance_prob) {
                x = x_new;
                accepted_moves++;
            }

            energy += wf.Energy(x);
        }
        dati[i] = energy / (nsteps / nblocks);
    }

    auto dati_elaborati = elaboraDati(nblocks, dati);
    writePsiSquaredToFile("psi_squared.dat", x, wf.PsiSquared(x));

    return std::make_pair(dati_elaborati[0][999], dati_elaborati[2][999]);
}

SimulatedAnnealingResult SimulatedAnnealing::run() {
    std::vector<double> dati(SA_steps);

    double sigma = rand.Rannyu(0, 1);
    double mu = rand.Rannyu(-1, 1);
    auto result = perform(sigma, mu);
    double energy_prev = result.first;
    double energy = 0;
    double error = 0;
    double current_temperature = initial_temperature;
    int t = 0;

    while (current_temperature > temp_min) {
        current_temperature = initial_temperature * pow(cooling_rate, t);
        std::string filename = "temp_" + std::to_string(current_temperature) + ".dat";
        initializeFile(filename);

        for (int step = 0; step < SA_steps; step++) {
            double candidate_sigma = sigma + rand.Rannyu(-0.1, 0.1);
            double candidate_mu = mu + rand.Rannyu(-0.1, 0.1);

            auto result = perform(candidate_sigma, candidate_mu);
            energy = result.first;
            error = result.second;

            writeDataToFile(filename, energy_prev, mu, sigma, error);

            double prob = exp(-(energy - energy_prev) / current_temperature);
            if (rand.Rannyu() < prob) {
                sigma = candidate_sigma;
                mu = candidate_mu;
                energy_prev = energy;
            }
        }
        t++;
    }

    SimulatedAnnealingResult final_result;
    final_result.sigma = sigma;
    final_result.mu = mu;
    final_result.last_energy = energy_prev;

    performMetropolisOptimized(sigma, mu);

    return final_result;
}

void SimulatedAnnealing::performMetropolisOptimized(double sigma, double mu) {
    int nsteps = 1000000;
    int nblocks = 10000;
    double step_size = 2.5;
    std::vector<double> dati(nblocks);
    double x = 0.0;
    int accepted_moves = 0;

    WaveFunction wf(sigma, mu);

    for (int i = 0; i < nblocks; i++) {
        double energy = 0;
        for (int j = 0; j < nsteps / nblocks; j++) {
            double x_new = x + rand.Rannyu(-step_size, step_size);
            double ratio = wf.PsiSquared(x_new) / wf.PsiSquared(x);
            double acceptance_prob = std::min(1.0, ratio);

            if (rand.Rannyu() < acceptance_prob) {
                x = x_new;
                accepted_moves++;
            }

            energy += wf.Energy(x);
        }
        dati[i] = energy / (nsteps / nblocks);
    }

    auto dati_elaborati = elaboraDati(nblocks, dati);

    std::ofstream outfile("optimized_energy.dat");
    if (outfile.is_open()) {
        outfile << std::left << std::setw(10) << "#ENERGY" << " "
                << std::setw(10) << "ERROR" << std::endl;
        for (int i = 0; i < nblocks; ++i) {
            outfile << std::left << std::setw(10) << dati_elaborati[0][i] << " "
                    << std::setw(10) << dati_elaborati[2][i] << std::endl;
        }
        outfile.close();
    } else {
        std::cerr << "Errore nell'apertura del file optimized_energy.dat" << std::endl;
    }
}

void SimulatedAnnealing::initializeFile(const std::string& filename) {
    std::ofstream outfile(filename);
    if (outfile.is_open()) {
        outfile << std::left << std::setw(10) << "#MU" << " "
                << std::setw(10) << "SIGMA" << " "
                << std::setw(10) << "ENERGY" << " "
                << std::setw(10) << "ERROR" << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void SimulatedAnnealing::writeDataToFile(const std::string& filename, double energy, double mu, double sigma, double error) {
    std::ofstream outfile(filename, std::ios::app);
    if (outfile.is_open()) {
        outfile << std::left << std::setw(10) << mu << " "
                << std::setw(10) << sigma << " "
                << std::setw(10) << energy << " "
                << std::setw(10) << error << std::endl;
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

std::vector<std::vector<double>> SimulatedAnnealing::elaboraDati(int numblocchi, const std::vector<double>& media) {
    std::vector<std::vector<double>> risultato(3, std::vector<double>(numblocchi));
    std::vector<double> sumprog(numblocchi);
    std::vector<double> vec_media2(numblocchi);

    for (int k = 0; k < numblocchi; k++) {
        vec_media2[k] = media[k] * media[k];
    }

    for (int i = 0; i < numblocchi; i++) {
        for (int j = 0; j <= i; j++) {
            sumprog[i] += media[j];
            risultato[1][i] += vec_media2[j];
        }
        risultato[0][i] = sumprog[i] / (i + 1);
        risultato[1][i] /= (i + 1);
        if (i > 0)
            risultato[2][i] = sqrt((risultato[1][i] - risultato[0][i] * risultato[0][i]) / i);
        else
            risultato[2][i] = 0;
    }

    return risultato;
}

void SimulatedAnnealing::initializePsiSquaredFile(const std::string& filename) {
    std::ofstream outfile(filename, std::ios::trunc);
    if (outfile.is_open()) {
        outfile << std::left << std::setw(10) << "#X" << " "
                << std::setw(10) << "PSI_SQUARED" << std::endl;
        outfile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}

void SimulatedAnnealing::writePsiSquaredToFile(const std::string& filename, double x, double psi_squared) {
    std::ofstream outfile(filename, std::ios::app);
    if (outfile.is_open()) {
        outfile << std::left << std::setw(10) << x << " "
                << std::setw(10) << psi_squared << std::endl;
        outfile.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}
