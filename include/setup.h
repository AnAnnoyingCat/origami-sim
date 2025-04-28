#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <json.hpp>
#include <fstream>
/**
 * @brief Set the up simulation parameters by reading them from the specified json file
 * 
 * @param dt 
 * @param vertexMass 
 * @param EA 
 * @param k_fold 
 * @param k_facet 
 * @param k_face 
 * @param zeta 
 */
void setup_simulation_params(std::string filename, double& dt, double& vertexMass, double& EA, double& k_fold, double& k_facet, double& k_face, double& zeta);

    /**
    * @brief 						    Sets up simulation with a simple square with one diagonal fold. Will support reading meshes in the future.
    * 
    * @param q 						    Generalized coordinates of vertex positions
    * @param qdot 					    Generalized velocities of vertex positions
    * @param x0 						Fixed point constraints
    * @param P 						    Fixed point constraint matrix
    * @param V 						    Vertices of simulation mesh
    * @param F 						    Faces of simulation mesh
    * @param alpha0 				    Nominal angles in flat state. Each row represents the angles of one face, in order α1_23, α2_31, α3_12.
    * @param E 						    Edges of simulation mesh
    * @param edge_target_angle 		    Target folded angle of each edge
    * @param l0 						Rest lengths of all springs
    * @param edge_adjacent_vertices     A vector of size 4 arrays storing for each fold the four relevant vertices in order RIGHT, LEFT, BEGIN, END.
    * @param k_axial                    A vector storing all axial stiffnesses
    * @param k_crease                   A vector storing all crease stiffnesses
    * @param EA                         Used in calculating axial stiffness
    * @param k_fold                     Fold stiffness constant
    * @param k_facet                    Facet stiffness constant
    * @param k_face                     Face stiffness constant
  */

 void setup_mesh(std::string filename, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &alpha0, Eigen::MatrixXi &E, Eigen::VectorXd& edge_target_angle, Eigen::VectorXd &l0, Eigen::MatrixXi& edge_adjacent_vertices, Eigen::VectorXd &k_axial, Eigen::VectorXd& k_crease, const double EA, const double k_fold, const double k_facet, const double k_face);
