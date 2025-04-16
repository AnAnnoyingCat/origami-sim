#include <setup.h>
#include <iostream>
#include <trig_helper_functions.h>

using json = nlohmann::json;

void setup_simulation_params(std::string filename, double& dt, double& vertexMass, double& EA, double& k_fold, double& k_facet, double& k_face, double& zeta){
	// Use nlohmann JSON to grab the simulation parameters from the specified filename.
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
		std::cout << "error: params file not found" << std::endl;     
	}
}

void setup(Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &x0, Eigen::SparseMatrix<double>& P, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &alpha0, Eigen::MatrixXi &E, Eigen::VectorXd& edge_target_angle, Eigen::VectorXd &l0, Eigen::VectorXd &k_axial, Eigen::VectorXd& k_crease, double& EA, double& k_fold, double& k_facet, std::vector<std::array<int, 4>>& edge_adjacent_vertices){
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
	edge_target_angle.resize(5);
	edge_target_angle << 0, 0, 0, 0, -90.0 * M_PI / 180;

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

void setup_mesh(std::string filename, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::VectorXd &x0, Eigen::SparseMatrix<double>& P, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &alpha0, Eigen::MatrixXi &E, Eigen::VectorXd& edge_target_angle, Eigen::VectorXd &l0, Eigen::MatrixXi edge_adjacent_vertices, Eigen::VectorXd &k_axial, Eigen::VectorXd& k_crease, const double EA, const double k_fold, const double k_facet, const double k_face){
	// TODO: import the mesh and fill all the relevant fields
	std::ifstream file(filename);

	if (file){
		json params = json::parse(file);

		// Assign Vertices
		if (params.contains("vertices_coords")){
			std::vector<std::vector<double>> coords = params["vertices_coords"].template get<std::vector<std::vector<double>>>();

			V.resize(coords.size(), 3);
			for (int i = 0; i < coords.size(); i++){
				V(i, 0) = coords[i][0];
				V(i, 1) = coords[i][1];
				V(i, 2) = 0;
			}
		} else {
			std::cout << "Missing parameter \"vertices_coords\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}

		// Assign faces
		if (params.contains("faces_vertices")){
			std::vector<std::vector<double>> faces = params["faces_vertices"].template get<std::vector<std::vector<double>>>();
			F.resize(faces.size(), 3);
			for (int i = 0; i < faces.size(); i++){
				F(i, 0) = faces[i][0];
				F(i, 1) = faces[i][1];
				F(i, 2) = faces[i][2];
			}
		} else {
			std::cout << "Missing parameter \"faces_vertices\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}
		
		// Assign resting face angles
		alpha0.resize(F.rows(), 3);
		for (int i = 0; i < F.rows(); i++){
			// Angle 1 of ith face is the angle at the first vertex, between the second and third
			getAngle(alpha0(i, 0), V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2)));
			// Angle 2 is at the second vertex, between the third and the first
			getAngle(alpha0(i, 1), V.row(F(i, 1)), V.row(F(i, 2)), V.row(F(i, 0)));
			// Angle 3 is at the third vertex between the first and the second
			getAngle(alpha0(i, 2), V.row(F(i, 2)), V.row(F(i, 0)), V.row(F(i, 1)));
		}

		// Assign edges
		if (params.contains("edges_vertices")){
			std::vector<std::vector<double>> edges = params["edges_vertices"].template get<std::vector<std::vector<double>>>();
			E.resize(edges.size(), 2);
			for (int i = 0; i < edges.size(); i++){
				E(i, 0) = edges[i][0];
				E(i, 1) = edges[i][1];
			}
		} else {
			std::cout << "Missing parameter \"edges_vertices\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}

		// Assign edge target angles
		if (params.contains("edges_foldAngle")){
			std::vector<double> angles = params["edges_foldAngle"].template get<std::vector<double>>();
			edge_target_angle.resize(angles.size());
			for (int i = 0; i < angles.size(); i++){
				edge_target_angle(i) = angles[i];
			}
		} else {
			std::cout << "Missing parameter \"edges_foldAngle\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}

		
		// Get edge type and calculate l0, kaxial and kcrease
		if (params.contains("edges_assignment")){
			std::vector<std::string> assignments = params["edges_assignment"].template get<std::vector<std::string>>();

			// For each edge calculate l0, k_axial and k_crease
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
				if (assignments[i] == "B"){
					k_crease(i) = -1;
				} else if (assignments[i] == "F"){
					k_crease(i) = k_facet * l0(i);
				} else { // Edge type == "M" or "V"
					k_crease(i) = k_fold * l0(i);
				}
			}
			
		} else {
			std::cout << "Missing parameter \"edges_assignment\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}
		

		// Initialize q and qdot
		q.resize(V.rows() * V.cols());
		qdot.resize(V.rows() * V.cols());

		Eigen::MatrixXd Vt = V.transpose();
		q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows() * Vt.cols());
		qdot.setZero();

		// Specify which vertices are fixed in space. fixedVerts(i) == 1 <=> Vertex i is fixed in space
		Eigen::VectorXd fixedVerts;
		fixedVerts.resize(q.size());
		// Only fix the first vertex in space so the simulation doesn't fly off or something
		fixedVerts.setZero();
		fixedVerts(0) = 1;

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

		// fill in edge_adjacent_vertices
		if (params.contains("edges_faces")){
			std::vector<std::vector<double>> edges_faces = params["edges_faces"].template get<std::vector<std::vector<double>>>();
			edge_adjacent_vertices.resize(E.size(), 4);
			for (int i = 0; i < E.rows(); i++){
				// Die erste Face ist glaub immer die rechte -> unique nummer der ersten face ist die rechte, unique nummer der zweiten face die linke, anfang, ende. sollte passen. dokumentation anschauen
			}
		} else {
			std::cout << "Missing parameter \"edges_faces\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}



	} else {
		std::cout << "There was an error reading your file!" << std::endl;
		if (filename == ""){
			std::cout << "No fold file specified, rendering default square" << std::endl;
		} else {
			std::cout << "fold file not found, rendering default square" << std::endl;
		}
		// If no file is specified, load default square
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
		edge_target_angle.resize(5);
		edge_target_angle << 0, 0, 0, 0, -90.0 * M_PI / 180;

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
		edge_adjacent_vertices.resize(5, 4);
		edge_adjacent_vertices.row(0) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(1) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(2) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(3) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(4) << 0, 3, 2, 1;
	}


	// Thing to set up::
	// Initial mesh (a square with diagonals)
    

}