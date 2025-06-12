#include "assemble_damping_stiffness.h"

void assemble_damping_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> k_Axial, double zeta) {
	int numEdges = E.rows();

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(numEdges * 12);
	int tripletcnt = 0;
	
	// Reuse H
	Eigen::Matrix<double, 6, 6> H;

	for (int currEdge = 0; currEdge < numEdges; currEdge++){
		double c = 2 * zeta * std::sqrt(k_Axial(currEdge));
		int v0 = E(currEdge, 0); 
		int v1 = E(currEdge, 1);

		Eigen::Vector3d qdot0 = qdot.segment<3>(3 * v0);
		Eigen::Vector3d qdot1 = qdot.segment<3>(3 * v1);
		
		// Top left block (just the diagonal)
		triplets[tripletcnt++] = {3 * v0, 3 * v0, -c};
		triplets[tripletcnt++] = {3 * v0 + 1, 3 * v0 + 1, -c};
		triplets[tripletcnt++] = {3 * v0 + 2, 3 * v0 + 2, -c};
		// Bottom right block (just the diagonal)
		triplets[tripletcnt++] = {3 * v1, 3 * v1, -c};
		triplets[tripletcnt++] = {3 * v1 + 1, 3 * v1 + 1, -c};
		triplets[tripletcnt++] = {3 * v1 + 2, 3 * v1 + 2, -c};
		// Top right block
		triplets[tripletcnt++] = {3 * v0, 3 * v1, c};
		triplets[tripletcnt++] = {3 * v0 + 1, 3 * v1 + 1, c};
		triplets[tripletcnt++] = {3 * v0 + 2, 3 * v1 + 2, c};
		// Bottom left block
		triplets[tripletcnt++] = {3 * v1, 3 * v0, c};
		triplets[tripletcnt++] = {3 * v1 + 1, 3 * v0 + 1, c};
		triplets[tripletcnt++] = {3 * v1 + 2, 3 * v0 + 2, c};
	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}