#include <Eigen/Dense>
#include <Eigen/Sparse>

/**
 * @brief Calculate kinetic energy of a single vertex
 * 
 * @param T 	Output: Kinetic energy
 * @param qdot 	Velocity of vertex
 * @param mass 	Mass of vertex
 */
void T_vert(double &T, Eigen::Ref<const Eigen::Vector3d> qdot, double mass);