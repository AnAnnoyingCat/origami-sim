#include <Eigen/Dense>
#include <trig_helper_functions.h>
#include <iostream>

/**
 * @brief Calculate the per fold force
 * 
 * @param f 		Resulting 12x1 Matrix with the forces on each point
 * @param q1 		Point 1 (Point to the RIGHT of edge)
 * @param q2 		Point 2 (Point to the LEFT of edge)
 * @param q3 		Point 3 (Origin point of the edge)
 * @param q4 		Point 4 (Point the edge points to)
 * @param k_crease	Stiffness for given crease
 * @param theta_Target The target theta of the given crease
 */
void F_crease(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1,  Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target, int creaseID);

/// @brief Function to set up the previous angle tracking fields
/// @param edge_target_angle 
void setup_prev_angle(int size);

void F_crease_maple(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target);