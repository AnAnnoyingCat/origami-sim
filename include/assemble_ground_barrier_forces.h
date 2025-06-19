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
void assemble_ground_barrier_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance);