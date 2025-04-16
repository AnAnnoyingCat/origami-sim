#include <setup.h>
#include <iostream>

using json = nlohmann::json;

void setup_simulation_params(std::string filename, double& dt, double& vertexMass, double& EA, double& k_fold, double& k_facet, double& k_face, double& zeta){
	// Use nlohmann JSON to grab the simulation parameters from the specified filename. If something isn't specified, set it to default value.
	std::ifstream file(filename);

	if (file){
		json params = json::parse(file);

		if (params.contains("dt")){
			dt = params["dt"].template get<double>();
			std::cout << "json lib works and it got dt = " << dt << std::endl;
		} else {
			dt = 0.001;
		}
		if (params.contains("vertexMass")){
			vertexMass = params["vertexMass"].template get<double>();
		} else {
			vertexMass = 1;
		}
		if (params.contains("EA")){
			EA = params["EA"].template get<double>();
		} else {
			EA = 5e4;
		}
		if (params.contains("k_fold")){
			k_fold = params["k_fold"].template get<double>();
		} else {
			k_fold = 1e3;
		}
		if (params.contains("k_facet")){
			k_facet = params["k_facet"].template get<double>();
		} else {
			k_facet = 1e3;
		}
		if (params.contains("k_face")){
			k_face = params["k_face"].template get<double>();
		} else {
			k_face = 2e2;
		}
		if (params.contains("zeta")){
			zeta = params["zeta"].template get<double>();
		} else {
			zeta = 0.25;
		}

	} else {
		if (filename == ""){
			std::cout << "No params file specified, using default params" << std::endl;
		} else{
			std::cout <<"Params file no found" << std::endl;
		}
		dt = 0.001;          
		vertexMass = 1;      
		EA = 1.0 * 5e4;      
		k_fold = 1e3;        
		k_facet = 1e3;       
		k_face = 2e2;        
		zeta = 0.25;         
	}
}

void setup(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &x0, Eigen::SparseMatrix<double>& P, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &alpha0, Eigen::MatrixXi &E, Eigen::VectorXd& edge_theta, Eigen::VectorXd &l0, Eigen::VectorXd &k_axial, Eigen::VectorXd& k_crease, double& EA, double& k_fold, double& k_facet, std::vector<std::array<int, 4>>& edge_adjacent_vertices){
	// Initial mesh (a square with diagonals)
    V.resize(4, 3);
    V << 0, 0, 0,
         0, 1, 0,
         1, 0, 0,
         1, 1, 0;
    
	// Two triangle Faces
    F.resize(2, 3);
    F << 0, 1, 2,
		 3, 2, 1;

	alpha0.resize(2, 3);
	alpha0 << 90 * M_PI / 180, 45 * M_PI / 180, 45 * M_PI / 180,
			  90 * M_PI / 180, 45 * M_PI / 180, 45 * M_PI / 180;
	
	// Square with one diagonal Edge
	E.resize(5, 2);
	E << 0, 1,
		 1, 3,
		 3, 2,
		 2, 0,
		 2, 1;
	
	// Diagonal should be folded to 60 degrees
	edge_theta.resize(5);
	edge_theta << 0, 0, 0, 0, -90.0 * M_PI / 180;

	//Type of each edge. "B" = Border, "M" = Mountain, "V" = Valley, "F" = Flat / Facet
	std::vector<std::string> edge_type;
	edge_type.resize(5);
	edge_type[0] = "B";
	edge_type[1] = "B";
	edge_type[2] = "B";
	edge_type[3] = "B";
	edge_type[4] = "V";

	

	// Basic Edge Lengths and stiffnesses
	l0.resize(E.rows());
	k_axial.resize(E.rows());
	k_crease.resize(E.rows());
	for (int i = 0; i < E.rows(); i++){
		int v0 = E(i, 0);
		int v1 = E(i, 1);
		Eigen::Vector3d difference = V.row(v0) - V.row(v1);

		// Set up l0
		l0(i) = difference.norm();

		// Set up axial stiffness dependent on l0
		k_axial(i) = EA / l0(i);

		// Set up crease stiffness dependent on l0
		if (edge_type[i] == "B"){
			k_crease(i) = -1;
		} else if (edge_type[i] == "F"){
			k_crease(i) = k_facet * l0(i);
		} else { // Edge type == "M" or "V"
			k_crease(i) = k_fold * l0(i);
		}
	}

	// Initialize q to vertex positions and qdot to zero
	q.resize(V.rows() * V.cols());
	qdot.resize(V.rows() * V.cols());

	Eigen::MatrixXd Vt = V.transpose();
	q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows() * Vt.cols());
	qdot.setZero();

	// Specify which vertices are fixed in space. fixedVerts(i) == 1 <=> Vertex i is fixed in space
	Eigen::VectorXd fixedVerts;
	fixedVerts.resize(4);
	fixedVerts << 0, 1, 1, 0;

	// Create the projection matrix P 
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.resize(3 * fixedVerts.size());

	int fixedCnt = 0;
	int tripletCnt = 0;
	for (int i = 0; i < q.size() / 3; i++){
		if (fixedVerts(i) == 0){
			// Current position is not fixed
			int base = 3 * i;
			triplets[tripletCnt++] = {base, base, 1.0};
			triplets[tripletCnt++] = {base + 1, base + 1, 1.0};
			triplets[tripletCnt++] = {base + 2, base + 2, 1.0};
			fixedCnt++;
		}
	}

	P.resize(q.size(), q.size());
	P.setFromTriplets(triplets.begin(), triplets.end());

	// Copied from the assignment a2
	x0 = q - P.transpose()*P*q;
	//correct q and qdot so they are the right size
    q = P*q;
    qdot = P*qdot;



	// Add the edge adjacent vertices, used in the calculation for the crease constraints
	edge_adjacent_vertices.resize(5);
	edge_adjacent_vertices[0] = {-1, -1, -1, -1};
	edge_adjacent_vertices[1] = {-1, -1, -1, -1};
	edge_adjacent_vertices[2] = {-1, -1, -1, -1};
	edge_adjacent_vertices[3] = {-1, -1, -1, -1};
	edge_adjacent_vertices[4] = {0, 3, 2, 1};
}