#include "LeximaxMaster.h"

void LeximaxMaster::defineModelConVar(bool print) {
	env = new GRBEnv();
	model = new GRBModel(*env);
	env->set(GRB_StringParam_LogFile, "Generated Log-files/LeximaxMaster.log");

	// Define 'MIN_M' variable, which we will maximize iteratively
	MIN_M = GRBVar();
	char name[13];
	sprintf_s(name, "MIN_M");
	MIN_M = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, name);

	// First add all the constraints needed
	// We will initialize the model by adding this constraint for all agents in M
	C = new GRBConstr[K->M.size()];
	int counter = 0;
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			C[counter] = model->addConstr(-MIN_M >= 0.0);
			counter++;
		}
	}

	C_Sum1 = GRBConstr();
	GRBLinExpr lin = 0.0;
	C_Sum1 = model->addConstr(lin == 1);

	// Now define the columns
	columns = std::vector<GRBColumn>();
	columns.resize(K->S.size());
	for (int s = 0; s < K->S.size(); s++) {
		counter = 0;
		for (int i = 0; i < K->I.n; i++) {
			if (K->M[i] == 1) {
				// It only makes sense to include the agents in 'M'
				if (K->S[s].x[i] == 1) {
					columns[s].addTerm(1, C[counter]);
				}
				counter++;
			}
		}
		// And for the constraint that forces the sum of the weights to be equal to one
		columns[s].addTerm(1, C_Sum1);

		K->block_solution(K->S[s]);
	}

	model->write("Generated Formulations/LeximaxMaster.lp");
	// Define variables
	w = std::vector<GRBVar>(K->S.size());
	for (int i = 0; i < K->S.size(); i++) {
		char name_w[13];
		sprintf_s(name_w, "w_%i", i);
		w[i] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, columns[i], name_w);
	}

	GRBLinExpr obj = MIN_M;
	model->setObjective(obj, GRB_MAXIMIZE);

	dual_C = std::vector<double>();
	dual_C_Sum1 = std::vector<double>();

	prob_agents = std::vector<double>();
	for (int i = 0; i < K->I.n; i++) {
		if (K->M[i] == 1) {
			prob_agents.push_back(0);
		}
	}

	if (print) {
		model->write("Generated Formulations/LeximaxMaster.lp");
	}
}

LeximaxMaster::LeximaxMaster(IPSolver* K_in, bool print) {
	K = new IPSolver(K_in, print);
	K->model->write("Generated Formulations/LeximaxMaster.lp");

	defineModelConVar(print);
}

LeximaxMaster::~LeximaxMaster() {
	delete[] C;
	delete K;
	delete model;
	delete env;
}

