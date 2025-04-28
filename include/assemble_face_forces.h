#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <face_constraints.h>

/**
 * @brief Assembles the forces that each face exerts on its vertices in trying to keep its original angle
 * 
 * @param f 		The force vector to be added to
 * @param q 		Generalized coordinates
 * @param F 		Faces of the CP
 * @param alpha0	Resting angles of the faces (size m x 3 where m is the number of faces)
 * @param k_crease 	
 */
void assemble_face_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::MatrixXd> alpha0, double k_face);