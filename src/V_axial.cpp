#include "V_axial.h"

void V_spring_particle_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness){
	double len = (q0 - q1).norm();
	// Avoid len close to zero
	if (len < 1e-10){
		V = 0;
		return;
	}
	V = -0.5 * stiffness * std::pow(len - l0, 2);
}