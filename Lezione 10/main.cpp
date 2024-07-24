#include "mpi.h"
#include "random.h"
#include "utility.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <vector>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int numCities = 110;
	Random rand;
	rand.initialize(rank); // Inizializza il generatore di numeri casuali con MPI rank
	vector<City> cities = createCitiesFromDataFile("cap_prov_ita.dat");

	GeneticAlgorithm ga(cities, 1000, 0.2, 0.9, rand);
	ga.initializePopulation();

	int gen = 1000; 
	int M = 10; // Frequenza di sincronizzazione

	double bestFitness = numeric_limits<double>::max();
	Individual bestIndividual{vector<int>(numCities)};

	for (int i = 0; i < gen; i++) {
			ga.evolve();
			Individual localBest = ga.getBestIndividual();


				if (localBest.fitness < bestFitness) {
						bestFitness = localBest.fitness;
						bestIndividual = localBest;
				}

				if (i % M == 0) {
						vector<int> localBestPath = bestIndividual.path;
						vector<int> globalBestPath(numCities);
						vector<int> allPaths(size * numCities);

						MPI_Gather(localBestPath.data(), numCities, MPI_INT, allPaths.data(), numCities, MPI_INT, 0, MPI_COMM_WORLD);

						if (rank == 0) {
								double globalBestFitness = bestFitness;
								vector<int> globalBestIndividualPath = localBestPath;

								for (int p = 1; p < size; p++) {
										vector<int> tmpPath(allPaths.begin() + p * numCities, allPaths.begin() + (p + 1) * numCities);
										Individual tmpIndividual(tmpPath);
										tmpIndividual.calculateFitness(cities);

										if (tmpIndividual.fitness < globalBestFitness) {
												globalBestFitness = tmpIndividual.fitness;
												globalBestIndividualPath = tmpPath;
										}
								}

								bestFitness = globalBestFitness;
								bestIndividual.path = globalBestIndividualPath;
								bestIndividual.calculateFitness(cities);
						}

						MPI_Bcast(bestIndividual.path.data(), numCities, MPI_INT, 0, MPI_COMM_WORLD);
				}
		}

		if (rank == 0) {
				cout << "Best path found: " << endl;
				for (int cityIndex : bestIndividual.path) {
						cout << cityIndex << " ";
				}
				cout << endl;
				cout << "Best path length: " << bestIndividual.fitness << endl;
		}
				saveBestIndividualToFile("italy.dat", bestIndividual, cities);

		MPI_Finalize();
		return 0;
}


