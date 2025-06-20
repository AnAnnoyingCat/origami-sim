#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

/**
 * @brief Calculates the gravity force acting on each particle
 * 
 * @param f 	Vector to add the forces onto
 * @param g 	Gravity force vector
 * @param mass 	Per vertex vertex mass
 */
void assemble_gravity_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::Vector3d> g, double mass);