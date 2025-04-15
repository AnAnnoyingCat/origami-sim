#include <trig_helper_functions.h>
#include <iostream>

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