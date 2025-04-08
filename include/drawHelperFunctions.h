#include <Eigen/Dense>
/**
 * @brief Function to update the global vertex positions V for the Viewer to display
 * 
 * @param V Vertices of simulation mesh to update Matrix<double, n, 3>
 * @param q Generalized coordinates of vertices Vector<double, 3*n>
 */
void updateV(Eigen::MatrixXd& V, Eigen::VectorXd& q);