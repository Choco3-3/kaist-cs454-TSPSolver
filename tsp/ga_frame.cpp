#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "ga_frame.hpp"

// defintions
#define RATE_SCALE 		100

// main engine
int ga_main(std::vector <city_node> node_list, int population_size, unsigned long long int total_evaluations, int local_threshold, float elite_rate, float mutation_rate, float mutation_size, int warn_threshold) {

	// initialise generation
	int generation = 0;

	// set up srand seed
	srand(std::chrono::system_clock::now().time_since_epoch().count());

	// int evaluation counter for early exit/convergence
	unsigned long long int evaluation_counter = 0;

	// generate possible genes from node_list size
	std::vector <int> main_genome;
	for (int i = 1; i <= node_list.size(); i++)
		main_genome.push_back(i);

	// generate initial population
	std::vector <Individual> population;
	population.reserve(population_size);
	for (int i = 0; i < population_size; i++) {
		std::vector <int> init_sol = create_genome(main_genome);
		population.push_back(Individual(init_sol, node_list));
	}

	// begin optimisation of initial population
	evaluation_counter += local_opt(population, node_list, local_threshold);

	// cancel out evaluation counter if total_evaluations is == 0
	if (total_evaluations == 0)
		evaluation_counter = 0;

	//std::cout << "Optimisation Done" << std::endl;
	//std::cout << "Generation: " << generation << std::endl;
	//print_population(population);	

	int termination = 0;
	int warning = 0;
	while (evaluation_counter <= total_evaluations && !termination) {
		// create new generation population
		population = new_gen(population, node_list, elite_rate, mutation_rate, mutation_size);
		// increment to next generation
		generation++;
		//std::cout << "Generation: " << generation << std::endl;
		//print_population(population);

		// get termination criteria by getting best result from the new_generation
		double best_fitness = population[0].fitness;

		// begin optimisation of initial population
		evaluation_counter += local_opt(population, node_list, local_threshold);

		// cancel out evaluation counter if total_evaluations is == 0
		if (total_evaluations == 0)
			evaluation_counter = 0;

		//terminate when a minimum has been reached both global or minimum (shown by non-improvement of 2 opt swap of best from new gen)
		if (best_fitness <= population[0].fitness) {
			//std::cout << "Termination Trigger" << std::endl;
			// terminate if warning reached, else set up warning
			if (warning == warn_threshold)
				termination = 1;
			else
				warning++;
		}
		else {
			warning = 0;
		}
	}

	// once termination criteria is met
	// sort population in order for best result
	std::sort(population.begin(), population.end(), less_than_key());

	// write to solution file
	std::ofstream outfile("solution.csv");

	//std::cout << "Best Route: " << std::endl;
	// output data
	for (int k = 0; k < population[0].chromosome.size(); k++) {
		outfile << population[0].chromosome[k] << std::endl;
	}
	
	//std::cout << "fitness: " << population[0].fitness << std::endl;
	// print best distance
	std::cout << population[0].fitness << std::endl;

	return 0;
}

// generate random solution by shuffling seed genome
std::vector <int> create_genome(std::vector <int> main_genome) {
	// initialise random number generator engine
	// source: https://www.techiedelight.com/shuffle-vector-cpp/
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(std::begin(main_genome), std::end(main_genome), std::default_random_engine(seed));

	return main_genome;
}

// optimise the population using 2-opt local search
// LK heuristic is a better local search heuristic but is complicated and unwieldy for the
// already implemented data structures. will try again if time remains
unsigned long local_opt(std::vector <Individual> &population, std::vector <city_node> node_list, int max_calc) {
	unsigned long calculations = 0;
	for (int i = 0; i < population.size(); i++) {
		calculations += population[i].local_2_opt(node_list, max_calc);
		//std::cout << "Population: " << i << " Done" << std::endl;
	}
	//std::cout << std::endl;
	return calculations;
}

// generate new generation by running selection, crossover and mutation with percentages for each
std::vector <Individual> new_gen(std::vector <Individual> old_pop, std::vector <city_node> node_list, float elite_rate, float mutation_rate, float mutation_size) {
	// random generator
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());

	// create new population vector
	std::vector <Individual> new_pop;

	// sort old population so best fitness is at start
	std::sort(old_pop.begin(), old_pop.end(), less_than_key());

	// select elite chromosomes and add to new population
	// count automatically floors as an int
	int elite_count = (elite_rate / RATE_SCALE) * old_pop.size(); 
	//std::cout << "elite: " << elite_count << std::endl;
	for (int i = 0; i < elite_count; i++)
		new_pop.push_back(old_pop[i]);

	// initiate roulette wheel probability allocation
	std::vector <double> selection_wheel;
	double offset = 0;
	unsigned long long int total_fitness = 0;
	for (int i = 0; i < old_pop.size(); i++)
		total_fitness += old_pop[i].fitness;
	for (int i = 0; i < old_pop.size(); i++) {
		double local_fitness = offset + (old_pop[i].fitness / total_fitness);
		selection_wheel.push_back(local_fitness);
		offset = local_fitness;
	}
	//std::cout << "offset: " << offset << std::endl;

	// initiate crossover with rest of population based on roulette wheel to select parents
	// after crossover is done mutate
	for (int i = elite_count; i < old_pop.size(); i++) {
		//std::cout << "Trigger: " << i << std::endl;
		// generates doubles between 0 and 1
		std::uniform_real_distribution <> real_distr(0, 1);

		// select two chromosomes to mate with each other
		std::vector <int> select_list;
		for (int j = 0; j < 2; j++) {
			double r = real_distr(generator);
			// invert selection wheel choice as lower fitness is better
			for (int k = 0; k < old_pop.size(); k++) {
				if (r < selection_wheel[k]) {
					select_list.push_back(old_pop.size() - k - 1);
					break;
				}
			}
		}

		/*
		std::cout << "size: " << select_list.size() << std::endl;
		for (int i = 0; i < select_list.size(); i++)
			std::cout << "Choice: " << select_list[i] << std::endl;
		*/

		// crossover the chosen chromosomes
		std::vector <int> new_cross_chromosome = cross_chromosome(old_pop[select_list[0]].chromosome, old_pop[select_list[1]].chromosome);
		// mutate some of the chromosomes except the elites depending on probability by mutation rate
		new_cross_chromosome = mutate_chromosome(new_cross_chromosome, mutation_rate, mutation_size);
		// push new chromosome to the new_populatiuon
		new_pop.push_back(Individual(new_cross_chromosome, node_list));
		//std::cout << "new_pop_size: " << new_pop[i].chromosome.size() << std::endl;
	}

	// sort new population so best fitness is at start
	std::sort(new_pop.begin(), new_pop.end(), less_than_key());

	return new_pop; 
}


// crossover two parent chromosomes into new chromosome
std::vector <int> cross_chromosome(std::vector <int> chromo_1, std::vector <int> chromo_2) {
	// random generator
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	// generates ints between 0 and size of chromosome_1
	std::uniform_int_distribution <int> int_distr(0, chromo_1.size());

	// select subset positions for chromo_1 to be added to chromo_2
	// pos_list[0] refers to start and pos_list[1] refers to end
	std::vector <int> pos_list;
	for (int i = 0; i < 2; i++)
		pos_list.push_back(int_distr(generator));

	// switch numbers around if one is bigger than the other
	if (pos_list[0] > pos_list[1]) {
		int temp = pos_list[0];
		pos_list[0] = pos_list[1];
		pos_list[1] = temp;
	}
	//std::cout << "start: " << pos_list[0] << " " << "end: " << pos_list[1] << std::endl << std::endl;
	
	// grab substring between those two positions
	std::vector <int> sub_chromosome((chromo_1.begin() + pos_list[0]), (chromo_1.begin() + pos_list[1]));

	// cross over chromosome 2 into chromosome 1 into a new return
	std::vector <int> new_crossed_chromosome;

	// push elements not in sub_chromosome into new_chromosome
	int sub_chromo_count = 0;
	int chromo_2_count = 0;
	for (int i = 0; i < chromo_1.size(); i++) {
		// if sub_chromo slice position, push in subchromosome
		if (i >= pos_list[0] && i < pos_list[1]) {
			new_crossed_chromosome.push_back(sub_chromosome[sub_chromo_count]);
			sub_chromo_count++;
		}
		// otherwise add in chromosome_2 gene that is not in sub_chromo
		else {
			while (chromo_2_count < chromo_2.size()) {
				// if we find a gene that is not in the sub_chromosome, push it and break
				if (!(std::find(sub_chromosome.begin(), sub_chromosome.end(), chromo_2[chromo_2_count]) != sub_chromosome.end())) {
					new_crossed_chromosome.push_back(chromo_2[chromo_2_count]);
					chromo_2_count++;
					break;
				}
				// if element is in array, move to the next one
				chromo_2_count++;
			}
		}
	}
	
	return new_crossed_chromosome;
}

// mutate the chromosome to ensure some randomness
// the mutation rate determines how big the mutation will be by determining the chance of a 
// random (between size of) subvector that ebing repositioned to another random location after inversion
// a mutation rate of 100 is not recommended as it does practically nothing
// try to stay below 75
std::vector <int> mutate_chromosome(std::vector <int> target_chromo, float mutation_rate, float mutation_size) {
	// random generator
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());

	// early exit if mutation will not go ahead
	std::uniform_real_distribution <> real_distr(0, 1);
	double mutation_chance = (mutation_rate / RATE_SCALE);
	if (mutation_chance < real_distr(generator))
		return target_chromo;

	// otherwise, continue with mutation
	// get mutation size
	int mutation_count = (mutation_size / RATE_SCALE) * target_chromo.size();
	// generates ints between 0 and size of target chromosome minus mutation vector size
	std::uniform_int_distribution <int> int_distr(0, (target_chromo.size() - mutation_count));

	// select subset positions for vector and repositioning
	// pos_list[0] refers to start of the repositioned vector and pos_list[1] refers to new start of the repositioned vector
	std::vector <int> pos_list;
	for (int i = 0; i < 2; i++)
		pos_list.push_back(int_distr(generator));

	// grab substring between those two positions
	std::vector <int> sub_chromosome((target_chromo.begin() + pos_list[0]), (target_chromo.begin() + pos_list[0] + mutation_count));

	// mutate chromosome into a new return
	std::vector <int> new_mutate_chromosome;

	// invert new vector and recombine
	// push elements not in sub_chromosome into new_chromosome
	int sub_chromo_count = 0;
	int chromo_count = 0;
	for (int i = 0; i < target_chromo.size(); i++) {
		// if sub_chromo slice position, push in subchromosome
		if (i >= pos_list[1] && i < (pos_list[1] + mutation_count)) {
			new_mutate_chromosome.push_back(sub_chromosome[(sub_chromosome.size() - sub_chromo_count) - 1]);
			sub_chromo_count++;
		}
		// otherwise add in chromosome_2 gene that is not in sub_chromo
		else {
			while (chromo_count < target_chromo.size()) {
				// if we find a gene that is not in the sub_chromosome, push it and break
				if (!(std::find(sub_chromosome.begin(), sub_chromosome.end(), target_chromo[chromo_count]) != sub_chromosome.end())) {
					new_mutate_chromosome.push_back(target_chromo[chromo_count]);
					chromo_count++;
					break;
				}
				// if element is in array, move to the next one
				chromo_count++;
			}
		}
	}

	return new_mutate_chromosome;
}

// Class functions
Individual::Individual(std::vector <int> chromosome_init, std::vector <city_node> node_list) {
	chromosome = chromosome_init;
	fitness = calculate_fitness(chromosome, node_list);
}

// calculate fitness of the chromosome using methods from node_list
double Individual::calculate_fitness(std::vector <int> target_chromo, std::vector <city_node> node_list) {
	double fitness = 0;

	// -1 as arrays start at 0 but nodes start at 1
	for (int i = 0; i < (target_chromo.size() - 1); i++)
		fitness += node_list[(target_chromo[i])-1].get_length(node_list[(target_chromo[i+1])-1]);

	// add fitness for travelling back to source city
	fitness += node_list[(target_chromo.back())-1].get_length(node_list[(target_chromo.front())-1]);

	return fitness;
}

// 2_opt local search optimisation method
unsigned long Individual::local_2_opt(std::vector <city_node> node_list, int threshold) {
	int improvement = 1;
	int size = chromosome.size();
	int fitness_calc = 0;

	// continue 2_opt while improvements are below the convergence threshold of chromosome.size()
	// continue 2_opt until fitness_calc threshhold reached
	while ((improvement <= chromosome.size()) && (fitness_calc <= threshold)) {
		// get initial fitness from object
		double best_fitness = fitness;

		for (int i = 0; (i < (size - 1)) && (fitness_calc < threshold); i++) {
			for (int k = i + 1; (k < size) && (fitness_calc < threshold); k++) {
				std::vector <int> new_chromosome = opt_2_swap(i, k);
				double new_fitness = calculate_fitness(new_chromosome, node_list);
				fitness_calc++;

				// cancel out fitness counter if threshold is == 0
				if (threshold == 0)
					fitness_calc = 0;

				// replace old chromosome and fitness if new_fitness better than old
				if (new_fitness < best_fitness) {
					chromosome = new_chromosome;
					fitness = new_fitness;
					best_fitness = fitness;
					improvement = 0;
					break;
				}
			}
			if (!improvement)
				break;
		}
		improvement++;
	} 
	return fitness_calc;
}

// opt_2_swap function
std::vector <int> Individual::opt_2_swap(int i, int k) {
	int size = chromosome.size();
	// make new chromosome for copy
	std::vector <int> new_chromosome;

	// 1. take route[0] to route[i-1] and add them in order to new_route
	for (int j = 0; j <= (i - 1); j++)
		new_chromosome.push_back(chromosome[j]);

	// 2. take route[i] to route[k] and add them in reverse order to new_route
	int rev_iter = 0;
	for (int j = i; j <= k; j++) {
		new_chromosome.push_back(chromosome[(k - rev_iter)]);
		rev_iter++;
	}

	// 3. take route[k+1] to end and add them in order to new_route
	for (int j = (k + 1); j < size; j++)
		new_chromosome.push_back(chromosome[j]);

	return new_chromosome;
}

// DEBUG FUNCTIONS
void print_population(std::vector <Individual> population) {
	for (int i = 0; i < population.size(); i++) {
		std::cout << "Population: " << i << std::endl;
		/*
		for (int k = 0; k < population[i].chromosome.size(); k++) {
			std::cout << population[i].chromosome[k] << " ";
		}
		*/
		std::cout << std::endl;
		std::cout << "fitness: " << population[i].fitness << std::endl << std::endl;
	}
}

