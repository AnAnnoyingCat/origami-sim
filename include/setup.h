#include <iostream>
#include <Eigen/Dense>

/**
 * @brief 		Sets up simulation with a simple square with one diagonal fold. Will support reading meshes in the future.
 * 
 * @param q 	Generalized coordinates of vertex positions
 * @param qdot 	Generalized velocities of vertex positions
 * @param x0 	Fixed point constraints
 * @param V 	Vertices of simulation mesh
 * @param F 	Faces of simulation mesh
 * @param E 	Edges of simulation mesh
 * @param l0 	Rest lengths of all springs
 */
void setup(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &x0, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &E, Eigen::VectorXd &l0);