#ifndef SIMULATEDANNEALING_H
#define SIMULATEDANNEALING_H

#include "WaveFunction.h"

struct SimulatedAnnealingResult {
    double sigma;
    double mu;
    double last_energy;
};

class SimulatedAnnealing {
public:
    SimulatedAnnealing(int SA_steps, double initial_temperature, double cooling_rate, double temp_min, Random& rand);
    SimulatedAnnealingResult run();
    void initializePsiSquaredFile(const std::string& filename);
    void writePsiSquaredToFile(const std::string& filename, double x, double psi_squared);

private:
    int SA_steps;
    double initial_temperature;
    double cooling_rate;
    double temp_min;
    Random& rand;

    std::pair<double, double> perform(double sigma, double mu);
    void performMetropolisOptimized(double sigma, double mu);
    void initializeFile(const std::string& filename);
    void writeDataToFile(const std::string& filename, double energy, double mu, double sigma, double error);
    std::vector<std::vector<double>> elaboraDati(int numblocchi, const std::vector<double>& media);
    
};

#endif // SIMULATEDANNEALING_H
