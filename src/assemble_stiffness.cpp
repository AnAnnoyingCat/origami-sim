#include "assemble_stiffness.h"

void assemble_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, Eigen::Ref<const Eigen::VectorXd> k){
	int numEdges = E.rows();
	
	K.resize(V.rows() * 3, V.rows() * 3);
	K.setZero();

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(numEdges * 36);
	int tripletcnt = 0;

	for (int i = 0; i < numEdges; i++){
		int v0 = E(i, 0); 
		int v1 = E(i, 1);

		Eigen::Vector3d q0 = q.segment<3>(3 * v0);
		Eigen::Vector3d q1 = q.segment<3>(3 * v1);

		Eigen::Matrix<double, 6, 6> H;
		d2V_axial_dq2(H, q0, q1, l0(i), k(i));
		H = -1 * H;

		// Iterate over all four stiffness blocks and fill in the triplets
		for (int b1 = 0; b1 < 2; b1++){
			for (int b2 = 0; b2 < 2; b2++){
				Eigen::Matrix3d block = H.block<3, 3>(3 * b1, 3 * b2);

				// Iterate through the entire block
				for (int row = 0; row < 3; row++){
					for (int col = 0; col < 3; col++){
						triplets[tripletcnt++] = {3 * E(i, b1) + row, 3 * E(i, b2) + col, block(row, col)};
					}
				}
			}
		}
	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}