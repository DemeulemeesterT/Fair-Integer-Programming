#include "IPSolver.h"


struct Ufl_parameters {
	int m; // Number of customers
	int n; // Number of locations

	// Minimum and maximum fixed costs to generate from
	int f_min;
	int f_max;

	// Minimum and maximum allocation costs to generate from
	int c_min;
	int c_max;

	// The name of the file to which the data should be exported
	std::string filename;

	// Seed that will be used
	unsigned seed;
	
};

#pragma once
class Ufl
{
public:
// Code for Uncapacitated Facility Location problem
	// Utility of an agent is a linear function of the negative cost of the facility to which she is assigned.
	// More specifically, we chose u_j = \max_{ij} c_{ij} - \sum_i c_{ij} x_{ij}
		// where j is agent
		// i is facility
		// c_ij is cost of agent i being assigned to facility j
		// x_ij is 1 is agent i is assigned to facility j
	
	// Note that the exact choice of the constant value in this expression will not affect the results	
		// As we will, for each agent, compare the utility to the worst utility she receives in any of the optimal solutions.

	// The model is the general model for uncapacitated facility location
		// See, e.g., Fischetti, Ljubic, Sinnl (2017), Management Science, Section 2.1 with an extra constraint for the utility.

	// We will store the model as
		// X[j] = utility of agent j (which was called u_j above)
		// Y[i] = 1 if facility i is opened
		// M[ij] = 1 if agent j is assigned (Matched) to facility i.



// FUNCTIONS
	inst generate_formulation(bool export_inst, bool print);
		// Will generate an 'inst' object which contains the formulation

	void generate_instance(Ufl_parameters U, bool print);
		// Will generate an instance of the uncapacitated facility location problem
		// Using the Bilde-Krarup instance generator that is included in UflLib
			// https://resources.mpi-inf.mpg.de/departments/d1/projects/benchmarks/UflLib/BildeKrarup.html

	void export_to_dat(inst I, bool print);



// GENERAL VARIABLES
	int n_locations;
		// The number of locations

	int m_customers;
		// The number of customers

	std::vector<double> f;
		// The vector with the costs of opening a facility on location i

	std::vector<std::vector<double>> c;
		// c[i,j] is the cost of allocating customer j to facility located in i

	double c_alloc_max;
		// Maximum of c_alloc, used to compute utility.

	std::string filename;


// CONSTRUCTOR AND DESTRUCTOR
	Ufl(Ufl_parameters Ufl_input, bool print);
	~Ufl();
};

