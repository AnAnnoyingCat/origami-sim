#include <Eigen/Dense>

/**
 * @brief Calculate the per fold force
 * 
 * @param f 		Resulting 12x1 Matrix with the forces on each point
 * @param q0 		Point 1 (Point to the RIGHT of edge)
 * @param q1 		Point 2 (Point to the LEFT of edge)
 * @param q2 		Point 3 (Origin point of the edge)
 * @param q3 		Point 4 (Point the edge points to)
 * @param stiffness crease stiffness
 */
void F_crease(Eigen::Ref<Eigen::Matrix<double, 12, 1>> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, double stiffness);
