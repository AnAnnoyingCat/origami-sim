#include "face_constraints.h"
#include <iostream>

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
	Eigen::Vector3d dq1, dq2, dq3, n;
	double currentAlpha, r31, r21;

	getNormal(n, q1, q2, q3);
	getAngle(currentAlpha, q1, q2, q3);

	//Squared edge len
	r21 = (q2 - q1).squaredNorm();
	r31 = (q3 - q1).squaredNorm();
	
	// Calculate the forces like they do in the paper
	dq2 = -n.cross(q2 - q1) / r21;
	dq3 = n.cross(q3 - q1) / r31;

	dq1 = -dq2 - dq3;

	f.segment<3>(0) = dq1;
	f.segment<3>(3) = dq2;
	f.segment<3>(6) = dq3;

	f *= -k_face * (currentAlpha - alpha0);
	//std::cout << "currentAlpha: " << currentAlpha * 180 / M_PI << ", alpha0: " << alpha0 * 180 / M_PI << std::endl;
}