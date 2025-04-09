#include "make_mass_matrix.h"

void make_mass_matrix(Eigen::SparseMatrix<double> &M, Eigen::Ref<const Eigen::VectorXd> q, double mass){

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(q.size());
	M.resize(q.size(), q.size());
	
	for (int i = 0; i < q.size(); i++){
		triplets[i] = {i, i, mass};
	}
	M.setFromTriplets(triplets.begin(), triplets.end());
}