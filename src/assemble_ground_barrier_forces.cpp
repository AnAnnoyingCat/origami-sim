#include "assemble_ground_barrier_forces.h"


void assemble_ground_barrier_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance) {
	Eigen::Vector3d scratchpad_f;
	for (int currVert = 0; currVert < q.size() / 3; currVert++){
		get_barrier_force_for_vertex(scratchpad_f, q.segment<3>(3 * currVert), min_barrier_distance);
		f.segment<3>(3 * currVert) -= scratchpad_f;
	}
}

void get_barrier_force_for_vertex(Eigen::Vector3d &f, const Eigen::Vector3d q0, double d){
	double z = q0(2);
	double force = 0.0;
	if (z > 0.0 && z < d){
		force = - 2.0 *(z - d) * std::log(z / d) - std::pow((z - d), 2) / z;
	}
	f << 0, 0, force;
}