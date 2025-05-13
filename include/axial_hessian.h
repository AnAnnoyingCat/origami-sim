#include <Eigen/Dense>

/**
 * @brief Calculate the per axis force stiffness
 * 
 * @param f 		Resulting 6x6 matrix 
 * @param q0 		Generalized coordinates of q0
 * @param q1 		Generalized coordinates of q1
 * @param l0 		Undeformed length of the (spring)edge
 * @param stiffness Stiffness constant for this (spring)edge
 */
void dF_axial(Eigen::Ref<Eigen::Matrix<double, 6, 6>> Hessian, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness);
