#include <iostream>
#include <cstring>
 #include <unistd.h>
#include "matrix.hpp"
#include "ga_frame.hpp"

// main tsp engine
int main(int argc, char *argv[]) {
	//int population_size, unsigned long long int total_evaluations, int local_threshold, float elite_rate, float mutation_rate, float mutation_size, int warn_threshold
	int pop_size = 10;
	int opt, local_threshold, warn_threshold = 0;
	float elite_rate, mutation_rate, mutation_size = 0;
	unsigned long long int total_evaluations = 0;

	// get filename first after checking for validity
	std::string filename;
	if (argc >= 2 && argc <= 16)
		filename = argv[1];
	else {
		std::cout << "Usage: ./tsp filename [-pflerst parameter]" << std::endl;
		return -1;
	}

	// generate city node_list
	//std::string filename = "rl11849_test.tsp.txt";
	//std::string filename = "att48.tsp.txt";
	std::vector <city_node> node_list = create_node_vector(filename);

	// get options
	while ((opt = getopt(argc, argv, "p:f:l:e:r:s:t:")) != -1) {
		switch (opt) {
			case 'p':
				pop_size = std::stoi(optarg);
				//std::cout << "pop_size: " << pop_size << std::endl;
				break;

			case 'f':
				total_evaluations = std::stoull(optarg);
				//std::cout << "total_evaluations: " << total_evaluations << std::endl;
				break;

			case 'l':
				local_threshold = std::stoi(optarg);
				//std::cout << "local_threshold: " << local_threshold << std::endl;
				break;

			case 'e':
				elite_rate = std::stof(optarg);
				//std::cout << "elite_rate: " << elite_rate << std::endl;
				break;

			case 'r':
				mutation_rate = std::stof(optarg);
				//std::cout << "mutation_rate" << mutation_rate << std::endl;
				break;

			case 's':
				mutation_size = std::stof(optarg);
				//std::cout << "mutation_size: " << mutation_size << std::endl;
				break;

			case 't':
				warn_threshold = std::stoi(optarg);
				//std::cout << "warn_threshold: " << warn_threshold << std::endl;
				break;

			case '?':
				return -1;
		}
	}

	ga_main(node_list, pop_size, total_evaluations, local_threshold, elite_rate, mutation_rate, mutation_size, warn_threshold); 

	return 0;
}

// ./tsp rl11849_test.tsp.txt -p 10 -f 10000000 -l 25000 -e 30 -r 50 -s 50 -t 10