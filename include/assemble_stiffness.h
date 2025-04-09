#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <d2V_axial_dq2.h>
/**
 * @brief 		Assemble the sparse stiffness matrix of the system
 * 
 * @param K 	Output: The 3nx3n sparse stiffness matrix which is the negative hessian of the potential energy function
 * @param q 	Generalized Coordinates
 * @param qdot 	Generalized velocity 
 * @param V 	The nx3 matrix of undeformed vertex positions
 * @param E 	The mx2 connectivity matrix
 * @param l0 	The undeformed length of each edge
 * @param k 	The Stiffness of each edge
 */
void assemble_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     Eigen::Ref<const Eigen::VectorXd> k);