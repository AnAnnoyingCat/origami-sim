#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <axial_constraints.h>

/**
 * @brief Calculate the size 3*n vector of forces acting from collision with the ground
 * 
 * @param f 
 * @param q 
 * @param min_barrier_distance 
 */
void assemble_ground_barrier_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance);

void get_barrier_stiffness_for_vertex(double& stiffness, const Eigen::Vector3d q, double d);