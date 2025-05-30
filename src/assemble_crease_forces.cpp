#include "assemble_crease_forces.h"
#include <iostream>
#include <fstream>

void assemble_crease_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> edge_adjacent_vertices, Eigen::Ref<const Eigen::VectorXd> k_crease, Eigen::Ref<const Eigen::VectorXd> edge_target_angle){
	// Pre-allocate force and reuse the same memory for performance
	Eigen::Matrix<double, 12, 1> force;
	for (int currCrease = 0; currCrease < edge_target_angle.size(); currCrease++){
		// Check if current edge is a border edge or an undriven crease
		if (std::isnan(edge_target_angle(currCrease))){
			continue;
		}

		// Get current edge information and then call F_crease
		// This will add the force from the creases "wanting to be folded".
		int v_right = edge_adjacent_vertices(currCrease, 0);
		int v_left = edge_adjacent_vertices(currCrease, 1);
		int v_begin = edge_adjacent_vertices(currCrease, 2);
		int v_end = edge_adjacent_vertices(currCrease, 3);

		Eigen::Vector3d q1 = q.segment<3>(3 * v_right);
		Eigen::Vector3d q2 = q.segment<3>(3 * v_left);
		Eigen::Vector3d q3 = q.segment<3>(3 * v_begin);
		Eigen::Vector3d q4 = q.segment<3>(3 * v_end);

		F_crease(force, q1, q2, q3, q4, k_crease(currCrease), edge_target_angle(currCrease), currCrease);
		f.segment<3>(3 * v_right) += force.segment<3>(0);
		f.segment<3>(3 * v_left) += force.segment<3>(3);
		f.segment<3>(3 * v_begin) += force.segment<3>(6);
		f.segment<3>(3 * v_end) += force.segment<3>(9);
	}
}

