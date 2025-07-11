#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <parameters.h>

/**
 * @brief Calculates positional updates based on current forces in the system
 * 
 * @param q 				generalized coordinates
 * @param qdot 				generalized velocities
 * @param dt 				time step size
 * @param force 			force(result, q, qdot, (first time calling this function this iteration?)) returns all the forces in the system
 * @param tmp_force			scratch space to store the total force in
 */
template <typename FORCE>
inline void forward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, FORCE &force, Eigen::VectorXd tmp_force, SimulationData& simulationData, SimulationParams& simulationParams){
	
	// Get the currnent deformed positions 
	get_deformed_positions(simulationData, simulationParams);

	// Detect active collisions
	simulationData.collisions.build(simulationData.collision_mesh, simulationData.deformed_vertices, simulationParams.min_barrier_distance);
	// Set up friction specific collisiosn
	simulationData.friction_collisions.build(simulationData.collision_mesh, simulationData.deformed_vertices, simulationData.collisions, simulationData.barrier_potential, simulationData.barrier_stiffness, simulationParams.mu);

	// Calculate all the forces
	force(tmp_force, q, qdot, true);

	// Perform forward Euler update. May be replaced by another scheme later.
	qdot += tmp_force * dt;
	q += qdot * dt;
}
