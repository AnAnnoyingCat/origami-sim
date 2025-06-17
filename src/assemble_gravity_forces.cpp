#include "assemble_gravity_forces.h"

void assemble_gravity_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::Vector3d> g, double mass) {
	for (int i = 0; i < f.size() / 3; i++){
		f.segment<3>(3 * i) += mass * g;
	}
}