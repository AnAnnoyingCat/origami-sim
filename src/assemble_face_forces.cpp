#include "assemble_face_forces.h"

void assemble_face_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::MatrixXi> alpha0, double k_crease) {
	// Pre-allocate force and reuse the same memory for performance
	Eigen::Matrix<double, 9, 1> force;
	// TODO: FACE FORCES
}