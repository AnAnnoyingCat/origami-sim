#include <iostream>
#include <Eigen/Dense>
#include <setup.h>

void setup(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &x0, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &E, Eigen::VectorXd edge_theta, Eigen::VectorXd &l0, Eigen::VectorXd &k, double EA, std::vector<std::array<int, 4>> edge_adjacent_vertices){
	// Initial mesh (a square with diagonals)
    V.resize(4, 3);
    V << 0, 0, 0,
         0, 1, 0,
         1, 0, 0,
         1, 1, 0;
    
	// Two triangle Faces
    F.resize(2, 3);
    F << 2, 1, 0,
		 2, 3, 1;
	
	// Square with one diagonal Edge
	E.resize(5, 2);
	E << 0, 1,
		 1, 3,
		 3, 2,
		 2, 0,
		 2, 1;
	
	edge_theta.resize(5);
	edge_theta << 0, 0, 0, 0, 60;
	
	// Basic Edge Lengths and stiffnesses
	l0.resize(E.rows());
	k.resize(E.rows());
	for (int i = 0; i < E.rows(); i++){
		int v0 = E(i, 0);
		int v1 = E(i, 1);
		Eigen::Vector3d difference = V.row(v0) - V.row(v1);
		l0(i) = difference.norm();
		k(i) = EA / l0(i);
	}

	// Initialize q to vertex positions and qdot to zero
	q.resize(V.rows() * V.cols());
	qdot.resize(V.rows() * V.cols());

	Eigen::MatrixXd Vt = V.transpose();
	q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows() * Vt.cols());
	qdot.setZero();

	// Specify which vertices are fixed in space
	x0.resize(2);
	x0 << 1, 2;

	// Add the edge adjacent vertices, used in the calculation for the crease constraints
	edge_adjacent_vertices.resize(4);
	edge_adjacent_vertices[0] = {-1, -1, -1, -1};
	edge_adjacent_vertices[1] = {-1, -1, -1, -1};
	edge_adjacent_vertices[2] = {-1, -1, -1, -1};
	edge_adjacent_vertices[3] = {0, 3, 2, 1};
}