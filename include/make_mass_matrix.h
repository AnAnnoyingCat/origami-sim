#include <Eigen/Dense>
#include <Eigen/Sparse>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  mass - the mass of each particle in the mass-spring system.
//Output:
//  M - sparse mass matrix for mass-spring system
/**
 * @brief Makes the diagonal matrix with just the mass on the diagonals
 * 
 * @param M 	Resulting sparse matrix 	
 * @param q 	Generalized vertex coordinates
 * @param mass 	Mass for each vertex 
 */
void make_mass_matrix(Eigen::SparseMatrix<double> &M, Eigen::Ref<const Eigen::VectorXd> q, double mass);