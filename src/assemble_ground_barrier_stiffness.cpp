#include "assemble_ground_barrier_stiffness.h"

void assemble_ground_barrier_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance, double k_barrier) {
	int numVerts = q.size() / 3;
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(numVerts);
	double scratchpad_stiffness;

	for (int currVert = 0; currVert < numVerts; currVert++){
		get_barrier_stiffness_for_vertex(scratchpad_stiffness, q.segment<3>(3*currVert), min_barrier_distance, k_barrier);
		triplets[currVert] = {3 * currVert + 2, 3 * currVert + 2, -1.0 * scratchpad_stiffness};
	}
	K.setFromTriplets(triplets.begin(), triplets.end());
}

void get_barrier_stiffness_for_vertex(double& stiffness, const Eigen::Vector3d q, double d, double k_barrier) {
	double z = q(2);
	double z2 = std::pow(z, 2.0);
	stiffness = 0.0;
	if (z > 0.0 && z < d){
		stiffness = -(2.0 * z2 * std::log(z / d) + 3.0*z2 - 2.0*d*z - std::pow(d, 2.0)) / z2;
		stiffness *= k_barrier;
	}
}