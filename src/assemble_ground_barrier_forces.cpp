#include "assemble_ground_barrier_forces.h"


void assemble_ground_barrier_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance, SimulationParams& simulationParams) {
	Eigen::Vector3d scratchpad_f;
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