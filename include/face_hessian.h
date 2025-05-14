 #include <Eigen/Dense>
#include <trig_helper_functions.h>
/**
 * @brief Calculate the per face force
 * 
 * @param f 		Resulting 9x9 stiffnesses
 * @param q1 		Point 1 (Give points in clockwise order)
 * @param q2 		Point 2 
 * @param q3 		Point 3
 * @param alpha0	Nominal angle in the flat states of angles α1_23, α2_31 and α3_12
 * @param k_face	Stiffness for given face
 */
void dF_face(Eigen::Ref<Eigen::Matrix<double, 9, 9>> f, Eigen::Ref<const Eigen::Vector3d> q1,  Eigen::Ref<const Eigen::Vector3d> q2, Eigen::Ref<const Eigen::Vector3d> q3, Eigen::Ref<const Eigen::Vector3d> alpha0, double k_face);
