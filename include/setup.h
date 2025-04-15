#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>


 /**
  * @brief 							Sets up simulation with a simple square with one diagonal fold. Will support reading meshes in the future.
  * 
  * @param q 						Generalized coordinates of vertex positions
  * @param qdot 					Generalized velocities of vertex positions
  * @param x0 						Fixed point constraints
  * @param P 						Fixed point constraint matrix
  * @param V 						Vertices of simulation mesh
  * @param F 						Faces of simulation mesh
  * @param alpha0 					Nominal angles in flat state. Each row represents the angles of one face, in order α1_23, α2_31, α3_12.
  * @param E 						Edges of simulation mesh
  * @param edge_theta 				Target folded angle of each edge
  * @param l0 						Rest lengths of all springs
  * @param k_axial 					Per axis stiffness constant calculated based on l0
  * @param k_crease 				Per crease stiffness constant calculated based on l0 of the crease axis
  * @param EA 						E is young's modulus and A is the cross-sectional area of the beam - Used in calculating k_axial
  * @param k_fold 					Constant double base stiffness in calculation for Mountain (M) and Valley (V) folds
  * @param k_facet 					Constant double base stiffness in calculation for Facet (F) folds
  * @param edge_adjacent_vertices 	A vector of size 4 arrays storing for each fold the four relevant vertices in order RIGHT, LEFT, BEGIN, END.
  */
void setup(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &x0, Eigen::SparseMatrix<double>& P, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &alpha0, Eigen::MatrixXi &E, Eigen::VectorXd& edge_theta, Eigen::VectorXd &l0, Eigen::VectorXd &k_axial, Eigen::VectorXd& k_crease, double& EA, double& k_fold, double& k_facet, std::vector<std::array<int, 4>>& edge_adjacent_vertices);