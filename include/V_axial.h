#include <Eigen/Dense>

/**
 * @brief Calculate potential energy of an edge. Used for graphing total PE in the system.
 * 
 * @param V 		Output potential energy of this (spring)edge
 * @param q0 		Generalized coordinates of first vertex
 * @param q1 		Generalized coordinates of second vertex
 * @param l0 		Undeformed length of the (spring)edge
 * @param stiffness Stiffness constant for this (spring)edge
 */
void V_spring_particle_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness);