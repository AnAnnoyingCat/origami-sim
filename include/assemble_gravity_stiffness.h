#include <Eigen/Dense>
#include <Eigen/Sparse>

/**
 * @brief 		This solely exist if I ever want to use a more complicated gravity formula. This does nothing
 * 
 * @param K 	If it did something, it would return it here.
 */
void assemble_gravity_stiffness(Eigen::SparseMatrix<double> &K);