#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <axial_constraints.h>

/**
 * @brief Calculate the size 3*n vector of forces acting on each vertex
 * 
 * @param f 		Resulting size 3*n vector
 * @param q 		Generalized coordinates
 * @param E 		mx2 connectivity matrix containing the edges
 * @param l0 		Undeformed length of every (spring)edge
 * @param k_Axial 	Size m vector of the stiffness of every (spring)edge
 */
void assemble_edge_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, Eigen::Ref<const Eigen::VectorXd> k_Axial);