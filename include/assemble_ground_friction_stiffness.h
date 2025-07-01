#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <axial_constraints.h>
#include <iostream>
#include <parameters.h>

/**
 * @brief Calculate the size 3*n x 3*n matrix of friction stiffnesses
 * 
 * @param f 
 * @param q 
 * @param min_barrier_distance 
 */
void assemble_ground_friction_stiffness(Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance);

/**
 * @brief Get the per vertex size 3x3 friction stiffness
 * 
 * @param K 
 * @param q0 
 * @param d 
 */
void get_friction_stiffness_for_vertex(Eigen::SparseMatrix<double> &K, const Eigen::Vector3d q0, double d);