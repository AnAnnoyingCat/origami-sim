#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

/**
 * @brief Calculates positional updates based on current forces in the system
 * 
 * @param q 				generalized coordinates
 * @param qdot 				generalized velocities
 * @param dt 				time step size
 * @param axial_force 		axial_force(result, q) returns the axial forces in the system
 * @param tmp_force_axial 	scratch space to store the axial force in (reserved for performance reasons)
 * @param tmp_force_crease  scratch space to store the crease force in
 * @param tmp_force_face	scratch space to store the face force in
 * @param tmp_force			scratch space to store the total force in
 */
template <typename AXIAL_FORCE>
inline void forward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, AXIAL_FORCE &axial_force, Eigen::VectorXd &tmp_force_axial, Eigen::VectorXd &tmp_force_crease, Eigen::VectorXd &tmp_force_face, Eigen::VectorXd tmp_force){
	// Begin by calculating all the forces
	axial_force(tmp_force_axial, q);

	// Crease forces

	// Face force

	tmp_force = axial_force;

	qdot += tmp_force * dt;
	q += qdot * dt;
}
