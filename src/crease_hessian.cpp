#include "crease_hessian.h"

void dF_crease(Eigen::Ref<Eigen::Matrix<double, 12, 12>> Hessian, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target) {
	Hessian.setZero();

	// Handwritten
	double q11 = q1(0);
	double q12 = q1(1);
	double q13 = q1(2);
	double q21 = q2(0);
	double q22 = q2(1);
	double q23 = q2(2);
	double q31 = q3(0);
	double q32 = q3(1);
	double q33 = q3(2);
	double q41 = q4(0);
	double q42 = q4(1);
	double q43 = q4(2);

	// Matlab generated code
	
	for (int i = 0; i < 12; i++){
		for (int j = 0; j < 12; j++){
			// set the result
		}
	}
}