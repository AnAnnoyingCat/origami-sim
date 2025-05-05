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
	double current_theta = std::atan2((n1.cross(n2)).dot(crease_dir), n1.dot(n2));
	if (std::abs(theta_target) > 175.0 * M_PI / 180.0){
		// once paper is close to fully folded force the sign to stay around there
		current_theta = std::copysign(current_theta, theta_target);
	} 
	//std::cout << "current angle: " << current_theta * 180 / M_PI << ", target angle: " << theta_target * 180 / M_PI << std::endl;

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
