#include "axial_hessian.h"

void dF_axial(Eigen::Ref<Eigen::Matrix<double, 6, 6>> Hessian, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
	Hessian.setZero();

	// Handwritten
	double q01 = q0(0);
	double q02 = q0(1);
	double q03 = q0(2);
	double q11 = q1(0);
	double q12 = q1(1);
	double q13 = q1(2);

	// ccode generated
	double H_simplified[6][6];
	H_simplified[0][0] = stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q02*q02)+l0*(q03*q03)+l0*(q12*q12)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q02*q12*2.0-l0*q03*q13*2.0);
	H_simplified[0][1] = -l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[0][2] = -l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[0][3] = -stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q02*q02)+l0*(q03*q03)+l0*(q12*q12)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q02*q12*2.0-l0*q03*q13*2.0);
	H_simplified[0][4] = l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[0][5] = l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[1][0] = -l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[1][1] = stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q03*q03)+l0*(q11*q11)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q03*q13*2.0);
	H_simplified[1][2] = -l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[1][3] = l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[1][4] = -stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q03*q03)+l0*(q11*q11)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q03*q13*2.0);
	H_simplified[1][5] = l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[2][0] = -l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[2][1] = -l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[2][2] = stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q02*q02)+l0*(q11*q11)+l0*(q12*q12)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q02*q12*2.0);
	H_simplified[2][3] = l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[2][4] = l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[2][5] = -stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q02*q02)+l0*(q11*q11)+l0*(q12*q12)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q02*q12*2.0);
	H_simplified[3][0] = -stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q02*q02)+l0*(q03*q03)+l0*(q12*q12)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q02*q12*2.0-l0*q03*q13*2.0);
	H_simplified[3][1] = l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[3][2] = l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[3][3] = stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q02*q02)+l0*(q03*q03)+l0*(q12*q12)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q02*q12*2.0-l0*q03*q13*2.0);
	H_simplified[3][4] = -l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[3][5] = -l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[4][0] = l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[4][1] = -stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q03*q03)+l0*(q11*q11)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q03*q13*2.0);
	H_simplified[4][2] = l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[4][3] = -l0*stiffness*(q01-q11)*(q02-q12)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[4][4] = stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q03*q03)+l0*(q11*q11)+l0*(q13*q13)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q03*q13*2.0);
	H_simplified[4][5] = -l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[5][0] = l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[5][1] = l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[5][2] = -stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q02*q02)+l0*(q11*q11)+l0*(q12*q12)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q02*q12*2.0);
	H_simplified[5][3] = -l0*stiffness*(q01-q11)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[5][4] = -l0*stiffness*(q02-q12)*(q03-q13)*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0);
	H_simplified[5][5] = stiffness*1.0/pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)*(l0*(q01*q01)+l0*(q02*q02)+l0*(q11*q11)+l0*(q12*q12)-pow(pow(q01-q11,2.0)+pow(q02-q12,2.0)+pow(q03-q13,2.0),3.0/2.0)-l0*q01*q11*2.0-l0*q02*q12*2.0);

	// Handwritten
	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){
			Hessian(i, j) = H_simplified[i][j];
		}
	}
}