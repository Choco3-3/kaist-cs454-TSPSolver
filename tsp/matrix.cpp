#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

#include "matrix.hpp"

// defintions
#define SKIP_LINE 			6
#define DIMENSION_LINE 		3
#define DIMENSION_DELIMIT 	2
#define NODE_LINE_TOKENS 	3
#define NODE_LINE_CITY 		0
#define NODE_LINE_XPOS 		1
#define NODE_LINE_YPOS 		2

// string definitions
const std::string EOF_LINE = "EOF";

// generate vector list of all nodes
std::vector <city_node> create_node_vector(const std::string &filename) {
	// use namespace std for less tediousness
	using namespace std;

	string line;
	ifstream datafile;
	int node_dimension = 0;

	// open datafile for reading
	datafile.open(filename);
	if (datafile.is_open()) {
		// skip lines
		for (int i = 0; i < SKIP_LINE; i++) {
			getline(datafile, line);
			// get dimensions using some evil pointer arithmetic magic
			if (i == DIMENSION_LINE)
				node_dimension = parse_dimension(line);
		}

		//create vector for created nodes
		std::vector <city_node> node_list;
		node_list.reserve(node_dimension);

		// read in Euclidian Data with space delimiter
		while (getline(datafile, line)) {
			// parse each line
			city_node *newNode;

			// early break out if EOF reached
			if ((newNode = parse_node(line)) == NULL)
				break;

			// add newNode object to end of list
			//std::cout << newNode->city_no << std::endl;
			node_list.push_back(*newNode);

			// destroy pointer
			delete newNode;
		}

		// close file
		datafile.close();
		return node_list;
	}
	else {
		cout << "file could not be opened" << endl;
	}

	// essentially return NULL
	return std::vector <city_node>();
}

// STATIC FUNCTIONS OTHER THAN MAIN USE STD:: FOR THE SAKE OF matrix.hpp
// tokenise function for parsing (USES ' ' ONLY)
static std::vector <std::string> tokenise(const std::string &line) {
	std::vector <std::string> tokens;
	std::stringstream line_stream(line);
	std::string parse;

	while (getline(line_stream, parse, ' '))
		tokens.push_back(parse);

	return tokens;
}

// parse dimension using tokenise function
static int parse_dimension(const std::string &line) {
	std::vector <std::string> tokens = tokenise(line);
	return stoi(tokens[DIMENSION_DELIMIT]);
}

// parse nodes using tokensie function, early return NULL pointer if line is EOF
static city_node *parse_node(const std::string &line) {
	std::vector <std::string> tokens = tokenise(line);

	//early check for EOF
	if (tokens.size() != NODE_LINE_TOKENS) {
		// check tokens not empty first to avoid potential segfault
		if (tokens.size() == 1 && tokens[0] == EOF_LINE) 
			return NULL;
		
		std::cout << "invalid number of tokens or no EOF detected" << std::endl;
		return NULL;
	}

	// parse node
	city_node *newNode = new city_node();
	newNode->city_no = stoi(tokens[NODE_LINE_CITY]);
	newNode->x_pos = stod(tokens[NODE_LINE_XPOS]);
	newNode->y_pos = stod(tokens[NODE_LINE_YPOS]);

	return newNode;
}