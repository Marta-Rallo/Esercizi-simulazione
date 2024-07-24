#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <numeric>
#include <unordered_set>
#include "random.h"

using namespace std;

// Define the City class
class City {
public:
		double x, y;

		City(double x, double y) : x(x), y(y) {}

		double distanceTo(const City& other_city) const {
				return sqrt(pow(x - other_city.x, 2) + pow(y - other_city.y, 2));
		}
};

// Define the TSP Individual class
class Individual {
public:
		vector<int> path;
		double fitness;

		Individual(const vector<int>& path) : path(path), fitness(0) {}

		void calculateFitness(const vector<City>& cities) {
				fitness = 0;
				for (size_t i = 0; i < path.size() - 1; i++) {
						fitness += cities[path[i]].distanceTo(cities[path[i + 1]]);
				}
				fitness += cities[path.back()].distanceTo(cities[path[0]]); // return to start
		}

		bool isValid(int numCities) {
				unordered_set<int> uniqueCities(path.begin(), path.end());
				return uniqueCities.size() == static_cast<std::unordered_set<int>::size_type>(numCities);
		}
};

// Define the Genetic Algorithm class
class GeneticAlgorithm {
public:
		vector<City> cities;
		vector<Individual> population;
		int populationSize;
		double mutationRate;
		double crossoverRate;
		Random rand;

		GeneticAlgorithm(const vector<City>& cities, size_t populationSize, double mutationRate, double crossoverRate)
				: cities(cities), populationSize(populationSize), mutationRate(mutationRate), crossoverRate(crossoverRate) {
				initializePopulation();
		}

		void initializePopulation() {
				vector<int> basePath(cities.size());
				iota(basePath.begin(), basePath.end(), 0); // initialize base path with 1, ..., N

				for (int i = 0; i < populationSize; i++) {
						random_shuffle(basePath.begin() + 1, basePath.end()); //c'è il +1 poichè voglio partire sempre dallo stesso punto per diminuire la degenerazione, random_shuffle mescola gli elementi del vettore
						Individual individual(basePath);
						individual.calculateFitness(cities);
						population.push_back(individual);
				}
		}

		Individual selectParent(double p) {
				// Sort the population based on fitness in descending order
				sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
						return a.fitness > b.fitness; // Assume higher fitness is better
				});

				// Generate a random number r in [0, 1)
				//double r = static_cast<double>(rand()) / RAND_MAX;
				double r = rand.Rannyu(0, 1);


				// Calculate the index j
				int M = population.size();
				int j = static_cast<int>(M * pow(r, p));

				// Return the selected individual
				return population[j];
		}

void shiftCities(Individual& individual, int start, int m, int n) {
		if (rand.Rannyu(0, 1) < mutationRate) {
				int N = individual.path.size();

				// Ensure valid range and not affecting the first city
				if (start < 1 || start >= N || m <= 0 || start + m > N) return;

				vector<int> temp(individual.path.begin() + start, individual.path.begin() + start + m);

				int shiftIndex = (start + n) % N;
				// Adjust shiftIndex to ensure it does not interfere with the start index
				if (shiftIndex < 1) shiftIndex += N;
				if (shiftIndex >= start && shiftIndex < start + m) return;

				// Perform the shift operation
				individual.path.erase(individual.path.begin() + start, individual.path.begin() + start + m);

				if (shiftIndex > start) {
						shiftIndex -= m; // Adjust shiftIndex to reflect the removal of elements before it
				}

				individual.path.insert(individual.path.begin() + shiftIndex, temp.begin(), temp.end());
		}
}


		void invertCities(Individual& individual, int start, int m) {
				if (rand.Rannyu(0, 1) < mutationRate) {
				int N = individual.path.size();
				if (start < 1 || start + m > N) return; // Ensure valid range and not affecting the first city

				reverse(individual.path.begin() + start, individual.path.begin() + start + m);
				}
		}

		void mutate(Individual& individual) {
				if (rand.Rannyu(0,1) < mutationRate) {
						int idx1 = static_cast<int>(rand.Rannyu(1,cities.size()));
						int idx2 = static_cast<int>(rand.Rannyu(1, cities.size()));

						swap(individual.path[idx1], individual.path[idx2]);
				}
		}

		void crossover(const Individual& parent1, const Individual& parent2, Individual& offspring1, Individual& offspring2) {
				if (rand.Rannyu(0, 1) < crossoverRate) {
						int cut = static_cast<int>(rand.Rannyu(1, parent1.path.size() - 1)); // Ensure valid cut point
						offspring1.path = vector<int>(parent1.path.begin(), parent1.path.begin() + cut);
						offspring2.path = vector<int>(parent2.path.begin(), parent2.path.begin() + cut);

						auto addRemainingCities = [](const vector<int>& parentPath, vector<int>& offspringPath) {
								unordered_set<int> citiesInOffspring(offspringPath.begin(), offspringPath.end());
								for (int city : parentPath) {
										if (citiesInOffspring.find(city) == citiesInOffspring.end()) {
												offspringPath.push_back(city);
												citiesInOffspring.insert(city);
										}
								}
						};

						addRemainingCities(parent2.path, offspring1.path);
						addRemainingCities(parent1.path, offspring2.path);
				} else {
						offspring1 = parent1;
						offspring2 = parent2;
				}
		}


void evolve() {
		vector<Individual> newPopulation;

		while (newPopulation.size() < static_cast<std::vector<Individual>::size_type>(populationSize)) {

				Individual parent1 = selectParent(0.2);
				Individual parent2 = selectParent(0.2);

				// Create offspring as copies of parents
				Individual offspring1(parent1.path);
				Individual offspring2(parent2.path);

				// Apply crossover operator
				crossover(parent1, parent2, offspring1, offspring2);

				// Apply mutation operator
				mutate(offspring1);
				mutate(offspring2);

				// Apply shiftCities operator
				int startShift = static_cast<int>(rand.Rannyu(1, offspring1.path.size() - 1)); // Random start position, excluding the first city
				int mShift = static_cast<int>(rand.Rannyu(1, offspring1.path.size() - startShift)); // Random number of cities to shift
				int nShift = static_cast<int>(rand.Rannyu(1, offspring1.path.size())); // Random shift amount

				shiftCities(offspring1, startShift, mShift, nShift);

				// Apply invertCities operator
				int startInvert = static_cast<int>(rand.Rannyu(1, offspring2.path.size())); // Random start position, excluding the first city
				int mInvert = static_cast<int>(rand.Rannyu(0, offspring2.path.size() - startInvert)); // Random number of cities to invert

				invertCities(offspring2, startInvert, mInvert);

				// Calculate fitness of offspring
				offspring1.calculateFitness(cities);
				offspring2.calculateFitness(cities);

				// Check validity of offspring
				if (offspring1.isValid(cities.size()) && offspring2.isValid(cities.size())) {
						newPopulation.push_back(offspring1);
						newPopulation.push_back(offspring2);
				}
		}

		// Replace the current population with the new population
		population = newPopulation;
}

 Individual getBestIndividual() {
				return *min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
						return a.fitness < b.fitness;
				});
		}
};

// Function to create cities positioned on a circumference
vector<City> createCitiesOnCircumference(int numCities, double radius, Random &rand) {
		vector<City> cities;
		for (int i = 0; i < numCities; ++i) {
				double angle = rand.Rannyu(0,2*M_PI);
				double x = radius * cos(angle);
				double y = radius * sin(angle);
				cities.push_back(City(x, y));
		}
		return cities;
}

vector<City> createCitiesInSquare(int numCities, double sideLength, Random& rand) {
		vector<City> cities;
		double halfSide = sideLength / 2.0;

		for (int i = 0; i < numCities; ++i) {
				double x = rand.Rannyu(-halfSide, halfSide);
				double y = rand.Rannyu(-halfSide, halfSide);
				cities.push_back(City(x, y));
		}

		return cities;
}

void saveBestIndividualToFile(const std::string& filename, const Individual& best, const std::vector<City>& cities) {
		std::ofstream outFile(filename);
		if (!outFile) {
				std::cerr << "Error opening file for writing: " << filename << std::endl;
				return;
		}

		for (int cityIndex : best.path) {
				outFile << cities[cityIndex].x << " " << cities[cityIndex].y << "\n";
		}

		outFile.close();
}


int main() {

		int numcities = 34;
		double radius = 1;
		Random rand;
		vector<City> cities = createCitiesOnCircumference(numcities, radius, rand);

		GeneticAlgorithm ga(cities, 1000, 0.2, 0.9);
		ga.initializePopulation();
		double bestFitness = numeric_limits<double>::max();
		Individual bestIndividual{vector<int>(numcities)};


		for (int generation = 0; generation < 1000; generation++) {
				ga.evolve();
				Individual best = ga.getBestIndividual();

				if (best.fitness <= bestFitness) {
						bestFitness = best.fitness;
						bestIndividual = best;
						cout << "Best path length = " << best.fitness << endl;
				}

		}
		saveBestIndividualToFile("best_individual_circle.dat", bestIndividual, cities);

		vector<City> cities_2 = createCitiesInSquare(numcities, 2, rand);
		GeneticAlgorithm ga_2(cities_2, 1000, 0.2, 0.9);
		ga_2.initializePopulation();
		double bestFitness_2 = numeric_limits<double>::max();
		Individual bestIndividual_2{vector<int>(numcities)};


		for (int generation = 0; generation < 1000; generation++) {
				ga_2.evolve();
				Individual best_2 = ga_2.getBestIndividual();

				if (best_2.fitness <= bestFitness_2) {
						bestFitness_2 = best_2.fitness;
						bestIndividual_2 = best_2;
						cout << "Best path length square = " << best_2.fitness << endl;
				}

		}
		saveBestIndividualToFile("best_individual_square.dat", bestIndividual_2, cities_2);

		return 0;
}
