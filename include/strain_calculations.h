#include <Eigen/Dense>

/// @brief      Helper function to calculate per face strain based on internal angles. Strain scales linearly from 0% to 100% angle deformation (e.g. angle of 20 deg deformed to angle of 40 deg)
/// @param C    Color matrix for each face
void calculateFaceAngleStrain(Eigen::MatrixXd& C, Eigen::MatrixXi& F, Eigen::VectorXd& q, Eigen::MatrixXd& alpha0);

/// @brief      Helper function to calculate per axis strain, visualized by recoloring the faces of the mesh. TODO this isn't implemented yet
/// @param C    Color matrix for each face
void calculateAxialDeformationStrain(Eigen::MatrixXd& C);