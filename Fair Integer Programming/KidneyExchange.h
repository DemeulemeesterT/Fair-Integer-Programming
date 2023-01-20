#include "IPSolver.h"


#pragma once
class KidneyExchange
{
public:
// FUNCTIONS
	void read_data(std::string name, bool print);

	void find_cycles(bool print);
		// Finds all cycles of length at most 3.

	inst generate_instance(bool export_inst, bool print);

	void export_to_dat(inst I, bool print);
		// Export the constructed instance to a .dat file


// GENERAL VARIABLES
	int n_pairs;
		// The number of agents

	int n_arcs;
		// The number of arcs in the graph

	std::vector<std::vector<bool>> A;
		// The arcs of the graph: A.[i][j] = 1 if (i,j) is in the graph
		// Arc (i,j) is in the graph if donor i is compatible with patient j

	std::vector<std::vector<bool>> C;
		// The set of all cycles of size at most 3 in the graph.
		// This is a matrix, 'C[c][i]' is 1 if agent i is included in cycle c

	std::string data_name;
		// Name of the instance


// CONSTRUCTOR AND DESTRUCTOR
	KidneyExchange(std::string name, bool print);
		// This constructor will create an instance based on the file 'name.INPUT'
		// This file describes the graph representing the possible exchange cycles
		// The format is inspired by the data used by:
			// Farnadi, G., St-Arnaud, W., Babaki, B., & Carvalho, M. (2021, May). Individual Fairness in Kidney Exchange Programs. 
				// In Proceedings of the AAAI Conference on Artificial Intelligence (Vol. 35, No. 13, pp. 11496-11505).

	~KidneyExchange();
};

