#include "assemble_friction_stiffness.h"

void assemble_friction_stiffness_IPC(Eigen::SparseMatrix<double> &K, SimulationParams &simulationParams, SimulationData &simulationData) {
    
	int full_size = simulationData.V.size() + simulationData.ground_V.size();

	// Fill velocities into a matrix
	Eigen::MatrixXd qdot_matrix_full(full_size, 3);
	qdot_matrix_full.setZero();
	for (int i = 0; i < simulationData.V.rows(); i++){
		qdot_matrix_full.row(i) = simulationData.qdot.segment<3>(3 * i);
	}

	// Calculate friction hessian
	Eigen::SparseMatrix<double> friction_potential_hess = simulationData.friction_potential.hessian(simulationData.friction_collisions, simulationData.collision_mesh, qdot_matrix_full);

	const int dim = simulationData.V.rows() * 3;

	// Extract top left block from the full hessian (we only apply forces to the model, not to the floor)
	std::vector<Eigen::Triplet<double>> triplets;
    for (int k = 0; k < friction_potential_hess.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(friction_potential_hess, k); it; ++it) {
            if (it.row() < dim && it.col() < dim) {
                triplets.emplace_back(it.row(), it.col(), it.value());
            }
        }
    }
	Eigen::SparseMatrix<double> hess_block(dim, dim);
    hess_block.setFromTriplets(triplets.begin(), triplets.end());

	K -= simulationParams.k_barrier * hess_block;
}