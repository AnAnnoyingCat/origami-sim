#include "make_mass_matrix.h"

void make_mass_matrix(Eigen::SparseMatrix<double> &M, Eigen::Ref<const Eigen::VectorXd> q, double mass){
	int numVert = q.size() / 3;
	
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(numVert);
	
	for (int i = 0; i < numVert; i++){
		triplets[i] = {i, i, mass};
	}
	M.setFromTriplets(triplets.begin(), triplets.end());
}