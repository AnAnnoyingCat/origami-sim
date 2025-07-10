#include "assemble_friction_stiffness.h"

void assemble_friction_stiffness_IPC(Eigen::SparseMatrix<double> &K, SimulationParams &simulationParams, SimulationData &simulationData) {
    // make a qdot_full with both model and ground vertices
	Eigen::VectorXd qdot_full(simulationData.qdot.size() + 3 * simulationData.ground_V.rows());
    
	qdot_full.setZero();
	qdot_full.head(simulationData.qdot.size()) = simulationData.qdot;
    std::cout << "=============== Stiffness ============" << std::endl;
    std::cout << "K.size(): " << K.size() << std::endl;
	std::cout << "deformed_vertices.rows(): " << simulationData.deformed_vertices.rows() << std::endl;
	std::cout << "collision_mesh.num_vertices(): " << simulationData.collision_mesh.num_vertices() << std::endl;
	std::cout << "qdot_full : " << qdot_full.transpose() << std::endl;

	Eigen::SparseMatrix<double> friction_potential_hess = simulationData.friction_potential.hessian(simulationData.friction_collisions, simulationData.collision_mesh, qdot_full);

	const int dim = simulationData.V.rows() * 3;

	// Extract top left block from the full hessian
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