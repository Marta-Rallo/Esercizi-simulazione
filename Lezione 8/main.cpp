#include <iostream>
#include "SimulatedAnnealing.h"

int main() {
    int SA_steps = 1500;
    double initial_temperature = 2;
    double cooling_rate = 0.9;
    double temp_min = 0.01;

    Random rand;
    SimulatedAnnealing sa(SA_steps, initial_temperature, cooling_rate, temp_min, rand);
    sa.initializePsiSquaredFile("psi_squared.dat");

    SimulatedAnnealingResult final_result = sa.run();

    std::cout << "Best parameters: sigma = " << final_result.sigma << ", mu = " << final_result.mu << std::endl;
    std::cout << "Best energy: " << final_result.last_energy << std::endl;

    return 0;
}
