#include "crease_constraints.h"
#include <iostream>

void F_crease(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target) {
	// Get all the relevant variables together
	double h1, h2;
	Eigen::Vector3d n1, n2;
	double alpha4_31, alpha3_14, alpha4_23, alpha3_42;
	
	geth(h1, q3, q4, q1);
	geth(h2, q3, q4, q2);
	getNormal(n1, q1, q4, q3);
	getNormal(n2, q2, q3, q4);
	getAngle(alpha4_31, q4, q3, q1);
	cot(alpha4_31);						// We only need the cot of the angles, so precompute it
	getAngle(alpha3_14, q3, q1, q4);
	cot(alpha3_14);
	getAngle(alpha4_23, q4, q2, q3);
	cot(alpha4_23);
	getAngle(alpha3_42, q3, q4, q2);
	cot(alpha3_42);

	// Calculate current fold angle theta 
    Eigen::Vector3d crease_dir = (q4 - q3).normalized();
    double cos_phi = n1.dot(n2);
    double sin_phi = (n1.cross(n2)).dot(crease_dir);
    double current_theta = M_PI - std::atan2(sin_phi, cos_phi);

	// Precompute some values
	Eigen::Vector3d n1h1 = n1 / h1;
	Eigen::Vector3d n2h2 = n2 / h2;
	
	double kdtheta = -k_crease * (current_theta - theta_target);

	// dθ/dp1
	f.segment<3>(0) = kdtheta * n1 / h1;

	// dθ/dp2
	f.segment<3>(3) = kdtheta * n2 / h2;
	
	// dθ/dp3
	f.segment<3>(6) = kdtheta * ((-alpha4_31 / (alpha3_14 + alpha4_31)) * n1h1 + (-alpha4_23 / (alpha3_42 + alpha4_23)) * n2h2);

	// dθ/dp4
	f.segment<3>(9) = kdtheta * ((-alpha3_14 / (alpha3_14 + alpha4_31)) * n1h1 + (-alpha3_42 / (alpha3_42 + alpha4_23)) * n2h2);\
}

/// Returns the angle at q0 between q1 (RIGHT) and q2 (LEFT) in alpha
void getAngle(double& alpha, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2){
	
	Eigen::Vector3d v1 = (q1 - q0).normalized();
	Eigen::Vector3d v2 = (q2 - q0).normalized();
	
	double cos_alpha = v1.dot(v2);

	alpha = std::acos(cos_alpha);
}

/// Returns the lever arm of the triangle perpendicular to the edge (q0, q1) in h
void geth(double& h, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2){
	Eigen::Vector3d base = q1 - q0;
	h = base.cross(q2 - q0).norm() / base.norm();
}

/// Returns the surface normal of triangle (q0, q1, q2) in n
void getNormal(Eigen::Vector3d& n, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2){
	n = ((q1 - q0).cross(q2 - q0)).normalized();
}

/// calculates the cotangent of alpha and stores it in alpha
void cot(double& alpha){
	// approximate it if alpha is very small
	if (std::fabs(alpha) < 1e-10){
		alpha = 1.0 / alpha;
	}

	double sin_alpha = std::sin(alpha);
    double cos_alpha = std::cos(alpha);

	if (std::fabs(sin_alpha) < 1e-10){
		std::cout << "WARNING: Very bad and small triangles" << std::endl;
	}

	alpha = cos_alpha / sin_alpha;

}