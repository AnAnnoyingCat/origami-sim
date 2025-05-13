#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <axial_hessian.cpp>
/**
 * @brief 		Assemble the sparse stiffness matrix of the system
 * 
 * @param K 	Output: The 3nx3n sparse stiffness matrix which is the negative hessian of the potential energy function
 * @param q 	Generalized Coordinates
 * @param V 	The nx3 matrix of undeformed vertex positions
 * @param E 	The mx2 connectivity matrix
 * @param l0 	The undeformed length of each edge
 * @param k 	The Stiffness of each edge
 */
void assemble_edge_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, Eigen::Ref<const Eigen::VectorXd> k);