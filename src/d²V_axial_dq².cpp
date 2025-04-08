#include <d²V_axial_dq².h>

void d2V_axial_dq2(Eigen::Ref<Eigen::Matrix<double, 6, 6>> H, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness){
	Eigen::Vector3d q = q1 - q0;
	double len = q.norm();

	// Avoid div by zero
	if (len < 1e-10){
		H.setZero();
		return;
	}
	Eigen::Matrix3d qqT = q * q.transpose();

	Eigen::Matrix3d dq0² = stiffness * (qqT / (std::pow(len, 2)) + (len - l0) * ((Eigen::Matrix3d::Identity() * (1.0 / len)) - (qqT / std::pow(len, 3))));
    
	H.block<3, 3>(0, 0) = dq0²;
	H.block<3, 3>(3, 0) = -dq0²;
	H.block<3, 3>(0, 3) = -dq0²;
	H.block<3, 3>(3, 3) = dq0²;
}