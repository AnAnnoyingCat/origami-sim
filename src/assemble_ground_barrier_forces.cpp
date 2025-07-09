#include "assemble_barrier_forces.h"



void assemble_ground_barrier_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance, SimulationParams& simulationParams, SimulationData simulationData) {
	Eigen::Vector3d scratchpad_f;
	Eigen::VectorXd f_old(f);
	for (int currVert = 0; currVert < q.size() / 3; currVert++){
		get_barrier_force_for_vertex(scratchpad_f, q.segment<3>(3 * currVert), min_barrier_distance, simulationParams.k_barrier);
		
		f.segment<3>(3 * currVert) -= scratchpad_f;
		
		if (currVert == 0){
			//std::cout << "Vertex 0 at height: " << q(2) << " experiencing force: " << scratchpad_f.transpose() * -1 << std::endl;
			if (q(2) < 0){
				std::cout << "Vertex 0 penetrated the ground, now at height: " << q(2) << std::endl;;
			}
		}
	}
	std::cout << "ground barrier forces: " << (f - f_old).transpose();
}

void get_barrier_force_for_vertex(Eigen::Vector3d &f, const Eigen::Vector3d q0, double d, double k_barrier){
	double z = q0(2);
	double force = 0.0;
	if (z > 0.0 && z < d){
		force = - 2.0 *(z - d) * std::log(z / d) - std::pow((z - d), 2) / z;
		force *= k_barrier;
	}
	f << 0, 0, force;
}

void assemble_barier_forces_IPC(Eigen::VectorXd &f, SimulationParams& simulationParams, SimulationData& simulationData, bool first_time){
	// Calculate the barrier potential derivative
	const ipc::BarrierPotential B(simulationParams.min_barrier_distance);
	Eigen::VectorXd barrier_potential_grad = B.gradient(simulationData.collisions, simulationData.collision_mesh, simulationData.deformed_vertices);

	if (simulationParams.enable_auto_k_barrier){
		// Update barrier stiffness
		if (!first_time){
			// Calculate current mindistance
			simulationData.curr_distance = simulationData.collisions.compute_minimum_distance(simulationData.collision_mesh, simulationData.deformed_vertices);
			
			// Calculate current barrier stiffness
			simulationData.barrier_stiffness = ipc::update_barrier_stiffness(
				simulationData.prev_distance, 
				simulationData.curr_distance, 
				simulationData.max_barrier_stiffness, 
				simulationData.barrier_stiffness, 
				simulationData.bbox_diagonal);
			
			// Update previous distance
			simulationData.prev_distance = simulationData.curr_distance;

		} else {
			// Get approximate scale of simulation
			simulationData.bbox_diagonal = ipc::world_bbox_diagonal_length(simulationData.deformed_vertices);

			// Make initial barrier stiffness guess
			simulationData.barrier_stiffness = ipc::initial_barrier_stiffness(
				simulationData.bbox_diagonal, 
				B.barrier(), 
				simulationParams.min_barrier_distance, 
				simulationParams.vertexMass, 
				f, 
				barrier_potential_grad, 
				simulationData.max_barrier_stiffness);
			
			// Calculate current minimum distance
			simulationData.prev_distance = simulationData.collisions.compute_minimum_distance(simulationData.collision_mesh, simulationData.deformed_vertices);
		}
		// Add the barrier stiffness * barrier force to the force vector
		f -= simulationData.barrier_stiffness * barrier_potential_grad.head(f.size());
	} else {
		f -= simulationParams.k_barrier * barrier_potential_grad.head(f.size());
	}
}