#include "T_vert.h"

void T_vert(double &T, Eigen::Ref<const Eigen::Vector3d> qdot, double mass){
	T = 0.5 * (qdot.transpose() * qdot).value();
}