#include "dV_axial_dq.h"

void dV_axial_dq(Eigen::Ref<Eigen::Matrix<double, 6, 1>> f, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness){
	Eigen::Vector3d q = q1 - q0;
	double len = q.norm();

	// Avoid div by zero
	if (len < 1e-10){
		f.setZero();
		return;
	}

	Eigen::Vector3d dVdq0 = stiffness * (len - l0) * q / len;

	// The two derivatives for q0 and q1 are exactly oppisite
	f.segment<3>(0) = dVdq0;
	f.segment<3>(3) = -dVdq0;
}