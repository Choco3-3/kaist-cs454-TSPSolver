// header guard
#ifndef GA_FRAME_H
#define GA_FRAME_H
// genetic algorithm framework header

#include <vector>
#include "matrix.hpp"

// class template for a possible solution (Individual) and required methods
// chromosome is the list of nodes chosen as the path
// fitness is the total distance of the trip
class Individual {
	public:
		// class/object attributes
		std::vector <int> chromosome;
		double fitness;

		// class object constructor
		Individual(std::vector <int> chromosome_init, std::vector <city_node> node_list);

		// class methods
		double calculate_fitness(std::vector <int> target_chromo, std::vector <city_node> node_list);
		unsigned long local_2_opt(std::vector <city_node> node_list, int threshold);
		std::vector <int> opt_2_swap(int i, int k);

};

// struct required for sort to work with object fitness attribute
// source: https://stackoverflow.com/questions/1380463/sorting-a-vector-of-custom-objects
struct less_than_key {
    inline bool operator() (const Individual &struct1, const Individual &struct2) {
        return (struct1.fitness < struct2.fitness);
    }
};

// main genetic algorithm engine
int ga_main(std::vector <city_node> node_list, int population_size, unsigned long long int total_evaluations, int local_threshold, float elite_rate, float mutation_rate, float mutation_size, int warn_threshold);

// generate random genome from size of node_list
std::vector <int> create_genome(std::vector <int> main_genome);

// optimise the population using 2-opt local search
unsigned long local_opt(std::vector <Individual> &population, std::vector <city_node> node_list, int max_calc);

// generate new generation by running selection, crossover and mutation
std::vector <Individual> new_gen(std::vector <Individual> old_pop, std::vector <city_node> node_list, float elite_rate, float mutation_rate, float mutation_size);

// crossover two parent chromosomes into new chromosome
std::vector <int> cross_chromosome(std::vector <int> chromo_1, std::vector <int> chromo_2);

// mutate the chromosome to ensure some randomness
std::vector <int> mutate_chromosome(std::vector <int> target_chromo, float mutation_rate, float mutation_size);

// debug functions
void print_population(std::vector <Individual> population);

#endif /* GA_FRAME_H */