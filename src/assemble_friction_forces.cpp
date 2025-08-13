#include "assemble_friction_forces.h"

void assemble_friction_forces_IPC(Eigen::VectorXd &f, SimulationParams &simulationParams, SimulationData &simulationData) {
	
	int full_size = simulationData.V.rows() + simulationData.ground_V.rows();

	// Fill velocities into a matrix
	Eigen::MatrixXd qdot_matrix_full(full_size, 3);
	qdot_matrix_full.setZero();
	for (int i = 0; i < simulationData.V.rows(); i++){
		qdot_matrix_full.row(i) = simulationData.qdot.segment<3>(3 * i);
	}

	// Calculate the friction gradient and subtract it from the forces
	Eigen::VectorXd friction_potential_grad = simulationData.friction_potential.gradient(simulationData.friction_collisions, simulationData.collision_mesh, qdot_matrix_full);
	f -= friction_potential_grad.head(f.size());
}