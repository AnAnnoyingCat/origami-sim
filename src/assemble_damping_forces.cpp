#include "assemble_damping_forces.h"

void assemble_damping_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> k_Axial, double zeta){
	// Pre allocate working memory
	Eigen::Vector3d force;
	for (int currEdge = 0; currEdge < E.rows(); currEdge++){
		double c = 2 * zeta * std::sqrt(k_Axial(currEdge));

		int v0 = E(currEdge, 0);
		int v1 = E(currEdge, 1);

		Eigen::Vector3d qdot0 = qdot.segment<3>(3 * v0);
		Eigen::Vector3d qdot1 = qdot.segment<3>(3 * v1);

		Eigen::Vector3d damping_force = c * (qdot1 - qdot0);
		f.segment<3>(3 * v0) += damping_force;
		f.segment<3>(3 * v1) -= damping_force;
	}
}