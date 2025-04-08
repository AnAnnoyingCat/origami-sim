#include "assemble_forces.h"

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, Eigen::Ref<const Eigen::VectorXd> k){
	// DONT forget to set f to zero first
	f.resize(q.size());
	f.setZero();

	for (int currEdge = 0; currEdge < E.rows(); currEdge++){
		// Get current edge information and call dV_axial_dq.
		// This will add the force from the potential energy of the edge deformation. Further forces may be added below.
		int v0 = E(currEdge, 0);
		int v1 = E(currEdge, 1);

		Eigen::Vector3d q0 = q.segment<3>(3 * v0);
		Eigen::Vector3d q1 = q.segment<3>(3 * v1);

		Eigen::Matrix<double, 6, 1> res;
		dV_axial_dq(res, q0, q1, l0(currEdge), k(currEdge));
		f.segment<3>(3 * v0) -= res.segment<3>(0);
		f.segment<3>(3 * v1) -= res.segment<3>(3);
	}
}