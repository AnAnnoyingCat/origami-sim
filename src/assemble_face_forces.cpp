#include "assemble_face_forces.h"
#include <iostream>

void assemble_face_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::MatrixXd> alpha0, double k_face) {
	// Pre-allocate force and reuse the same memory for performance
	Eigen::Matrix<double, 9, 1> force;
	for (int currentFace = 0; currentFace < F.rows(); currentFace++){
		int v1, v2, v3;
		Eigen::Vector3d q1, q2, q3, currAlpha0;

		v1 = F.row(currentFace)(0);
		v2 = F.row(currentFace)(1);
		v3 = F.row(currentFace)(2);
		
		q1 = q.segment<3>(3 * v1);
		q2 = q.segment<3>(3 * v2);
		q3 = q.segment<3>(3 * v3);

		currAlpha0 = alpha0.row(currentFace);

		F_face(force, q1, q2, q3, currAlpha0, k_face);
		
		f.segment<3>(3 * v1) += force.segment<3>(0);
		f.segment<3>(3 * v2) += force.segment<3>(1);
		f.segment<3>(3 * v3) += force.segment<3>(2);
	}
}