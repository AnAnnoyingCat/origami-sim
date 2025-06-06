#include <Eigen/Dense>
#include <trig_helper_functions.h>
#include <crease_constraints.h>

/**
 * @brief Calculate the per fold stiffness 
 * 
 * @param f 		Resulting 12x1 Matrix with the forces on each point
 * @param q1 		Point 1 (Point to the RIGHT of edge)
 * @param q2 		Point 2 (Point to the LEFT of edge)
 * @param q3 		Point 3 (Origin point of the edge)
 * @param q4 		Point 4 (Point the edge points to)
 * @param k_crease	Stiffness for given crease
 * @param theta_Target The target theta of the given crease
 */
void dF_crease(Eigen::Ref<Eigen::Matrix<double, 12, 12>> Hessian, Eigen::Ref<const Eigen::Vector3d> q1,  Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> q4, double k_crease, double theta_target);
