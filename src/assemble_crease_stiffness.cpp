#include "assemble_crease_stiffness.h"

void assemble_crease_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> edge_adjacent_vertices, Eigen::Ref<const Eigen::VectorXd> k_crease, Eigen::Ref<const Eigen::VectorXd> theta_target) {
	
	int numCreases = edge_adjacent_vertices.rows();

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(numCreases * 144);
	int tripletcnt = 0;
	
	// Reuse H
	Eigen::Matrix<double, 12, 12> H;

	for (int currCrease = 0; currCrease < numCreases; currCrease++){
		int v0 = edge_adjacent_vertices(currCrease, 0); 
		int v1 = edge_adjacent_vertices(currCrease, 1); 
		int v2 = edge_adjacent_vertices(currCrease, 2); 
		int v3 = edge_adjacent_vertices(currCrease, 3); 

		Eigen::Vector3d q1 = q.segment<3>(3 * v0);
		Eigen::Vector3d q2 = q.segment<3>(3 * v1);
		Eigen::Vector3d q3 = q.segment<3>(3 * v2);
		Eigen::Vector3d q4 = q.segment<3>(3 * v3);
		
		dF_crease(H, q1, q2, q3, q4, k_crease(currCrease), theta_target(currCrease));

		// Iterate over all 16 stiffness blocks and fill in the triplets
		for (int b1 = 0; b1 < 4; b1++){
			for (int b2 = 0; b2 < 4; b2++){
				Eigen::Matrix3d block = H.block<3, 3>(3 * b1, 3 * b2);

				// Iterate through the entire block
				for (int row = 0; row < 3; row++){
					for (int col = 0; col < 3; col++){
						triplets[tripletcnt++] = {3 * edge_adjacent_vertices(currCrease, b1) + row, 3 * edge_adjacent_vertices(currCrease, b2) + col, block(row, col)};
					}
				}
			}
		}
	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}

