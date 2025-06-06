#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <crease_hessian.h>

/**
 * @brief Calculate the size 3*n vector of forces acting on each vertex
 * 
 * @param K 							Stiffness matrix to be added onto
 * @param q 							Generalized coordinates
 * @param edge_adjacent_Vertices 		the size m vector of arrays of size 4 containing the triangles relevant for the current edge
 * @param k_Crease 						Size m vector of the stiffness of every crease, is -1 for border edges
 * @param curr_theta					Size m vector containing the current target theta for each crease
 */
void assemble_crease_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> edge_adjacent_vertices, Eigen::Ref<const Eigen::VectorXd> k_crease, Eigen::Ref<const Eigen::VectorXd> theta_target);