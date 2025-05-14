#include <Eigen/Dense>
#include <Eigen/Sparse>

/**
 * @brief Calculates the damping force of each vertex based on the velocity of its neighbours
 * 
 * @param f 		Vector where the force is returned into
 * @param qdot 		Generalized velocity of the system
 * @param E 		Edge connectivity matrix
 * @param k_Axial 	Per edge axial stiffness
 * @param zeta 		Parameter in the viscous damping coefficient
 */
void assemble_damping_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> k_Axial, double zeta);