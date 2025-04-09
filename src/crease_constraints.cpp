#include "crease_constraints.h"

void F_crease(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, double stiffness) {

}

/// Returns the angle at q0 between q1 (RIGHT) and q2 (LEFT) in alpha
void get_angle(double& alpha, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2){
	
	Eigen::Vector3d v1 = q1 - q0;
	Eigen::Vector3d v2 = q2 - q0;
	
	double cos_alpha = v1.dot(v2);

	alpha = std::acos(cos_alpha);
}

/// Returns the lever arm of the triangle perpendicular to the edge (q0, q1) in h
void geth(Eigen::Vector3d& h, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2){

}

/// Returns the surface normal of triangle (q0, q1, q2)
void getNormal(Eigen::Vector3d& h, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2){

}