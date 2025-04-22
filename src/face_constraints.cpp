#include "face_constraints.h"

void F_face(Eigen::Ref<Eigen::Matrix<double, 9, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> alpha0, double k_face) {
	// Each interior angle exerts forces on all three vertices of the triangle. 
	Eigen::Matrix<double, 9, 1> Ffacealpha1_23, Ffacealpha2_31, Ffacealpha3_12;
	
	calculateFfacealpha(Ffacealpha1_23, q1, q2, q3, alpha0[0], k_face);
	calculateFfacealpha(Ffacealpha2_31, q2, q3, q1, alpha0[1], k_face);
	calculateFfacealpha(Ffacealpha3_12, q3, q1, q2, alpha0[2], k_face);

	f = Ffacealpha1_23 + Ffacealpha2_31 + Ffacealpha3_12;
}

/// @brief  Assume alpha0 is alpha0q1_q2q3. 
void calculateFfacealpha(Eigen::Ref<Eigen::Matrix<double, 9, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, double alpha0, double k_face){
	Eigen::Vector3d dp1, dp2, dp3, n;
	double currentAlpha, squaredNorm31, squaredNorm21;

	getNormal(n, q1, q2, q3);
	getAngle(currentAlpha, q1, q2, q3);
	squaredNorm31 = (q3 - q1).squaredNorm();
	squaredNorm21 = (q2 - q1).squaredNorm();

	dp1 = (n.cross(q3 - q1)) / squaredNorm31;
	dp2 = (-n.cross(q3 - q1)) / squaredNorm31 + (n.cross(q2 - q1)) / squaredNorm21;
	dp3 = (-n.cross(q2 - q1)) / squaredNorm21;

	f.segment<3>(0) = dp1;
	f.segment<3>(3) = dp2;
	f.segment<3>(6) = dp3;

	f *= -k_face * (currentAlpha - alpha0);
}