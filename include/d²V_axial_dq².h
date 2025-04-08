#include <Eigen/Dense>

/**
 * @brief Calculates the 6x6 dense hessian matrix of potential energy of a (spring)edge
 * 
 * @param H 		Resulting 6x6 Matrix
 * @param q0 		Generalized coordinates of q0
 * @param q1 		Generalized coordinates of q1
 * @param l0 		Undeformed length of the (spring)edge
 * @param stiffness Stiffness constant for this (spring)edge
 */
void d2V_axial_dq2(Eigen::Ref<Eigen::Matrix<double, 6, 6>> H, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness);
