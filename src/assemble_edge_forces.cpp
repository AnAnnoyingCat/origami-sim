#include "assemble_edge_forces.h"

void assemble_edge_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, Eigen::Ref<const Eigen::VectorXd> k){

	for (int currEdge = 0; currEdge < E.rows(); currEdge++){
		// Get current edge information and then call dV_axial_dq.
		// This will add the force from the potential energy of the edge deformation.
		int v0 = E(currEdge, 0);
		int v1 = E(currEdge, 1);

		Eigen::Vector3d q0 = q.segment<3>(3 * v0);
		Eigen::Vector3d q1 = q.segment<3>(3 * v1);

		Eigen::Matrix<double, 6, 1> force;
		F_axial(force, q0, q1, l0(currEdge), k(currEdge));
		f.segment<3>(3 * v0) += force.segment<3>(0);
		f.segment<3>(3 * v1) += force.segment<3>(3);
	}
}