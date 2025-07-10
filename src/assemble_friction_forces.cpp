#include "assemble_friction_forces.h"

void assemble_friction_forces_IPC(Eigen::VectorXd &f, SimulationParams &simulationParams, SimulationData &simulationData) {
	

	// make a qdot_full with both model and ground vertices
	Eigen::VectorXd qdot_full(simulationData.qdot.size() + 3 * simulationData.ground_V.rows());
	qdot_full.setZero();
	qdot_full.head(simulationData.qdot.size()) = simulationData.qdot;
	std::cout << "=============== Friction ============" << std::endl;
	std::cout << "f.size(): " << f.size() << std::endl;
	std::cout << "Expected: " << 3 * simulationData.rest_positions.rows() << std::endl;
	std::cout << "deformed_vertices.rows(): " << simulationData.deformed_vertices.rows() << std::endl;
	std::cout << "collision_mesh.num_vertices(): " << simulationData.collision_mesh.num_vertices() << std::endl;
	std::cout << "qdot_full : " << qdot_full.transpose() << std::endl;

	std::cout << "qdot_full.size(): " << qdot_full.size() << std::endl;
	std::cout << "collision_mesh.num_vertices(): " << simulationData.collision_mesh.num_vertices() << std::endl;
	std::cout << "Number of collisions: " << simulationData.friction_collisions.size() << std::endl;

	Eigen::VectorXd friction_potential_grad = simulationData.friction_potential.gradient(simulationData.friction_collisions, simulationData.collision_mesh, qdot_full);
	
	std::cout << "friction_potential_grad.size(): " << friction_potential_grad.size() << std::endl;
	std::cout << "f.size(): " << f.size() << std::endl;
	f -= friction_potential_grad.head(f.size());
}