// header guard
#ifndef MATRIX_H
#define MATRIX_H
// matrix calculation functions header
#include <vector>
#include <cmath>
#include <string>

//class template for parsed node data objects
class city_node {
	public:
		int city_no;
		double x_pos;
		double y_pos;

		double get_length(const city_node &dest) {
			double dx = dest.x_pos - x_pos;
			double dy = dest.y_pos -y_pos;
			return std::sqrt((dx*dx) + (dy*dy)); 
		}
};

// node vector creation function
std::vector <city_node> create_node_vector(const std::string &filename);

// tokeniser for parsiong
static std::vector <std::string> tokenise(const std::string &line);
// dimension line parser
static int parse_dimension(const std::string &line);
// node line parser
static city_node *parse_node(const std::string &line);

#endif /* MATRIX_H */