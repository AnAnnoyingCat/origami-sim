#include <setup.h>
#include <iostream>
#include <trig_helper_functions.h>

using json = nlohmann::json;

void setup_simulation_params(std::string filename, double& dt, double& vertexMass, double& EA, double& k_fold, double& k_facet, double& k_face, double& zeta, bool& ENABLE_STRAIN_VISUALIZATION){
	// Use nlohmann JSON to grab the simulation parameters from the specified filename.
	std::ifstream file(filename);

	if (file){
		json params = json::parse(file);

		if (params.contains("dt")){
			dt = params["dt"].template get<double>();
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
		if (params.contains("visualize_strain")){
			ENABLE_STRAIN_VISUALIZATION = params["visualize_strain"].template get<bool>();
		} else {
			ENABLE_STRAIN_VISUALIZATION = false;
		}
		
	} else {
		std::cout << "error: params file not found" << std::endl;     
	}
}

void setup_mesh(std::string filename, Eigen::VectorXd &q, Eigen::VectorXd &qdot, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &alpha0, Eigen::MatrixXi &E, Eigen::VectorXd& edge_target_angle, Eigen::VectorXd &l0, Eigen::MatrixXi& edge_adjacent_vertices, Eigen::VectorXd &k_axial, Eigen::VectorXd& k_crease, const double EA, const double k_fold, const double k_facet, const double k_face){
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
				if (coords[i].size() == 3){
					V(i, 2) = coords[i][2];
				} else {
					V(i, 2) = 0;
				}
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
			std::cout << "we do have vertices coords. proof: " << V << std::endl;
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
			const auto& angles_json = params["edges_foldAngle"];
			edge_target_angle.resize(angles_json.size());

			for (int i = 0; i < angles_json.size(); i++){
				if (angles_json[i].is_null()){
					// null angle usually means border edge or facet crease. in any case, set 0.
					edge_target_angle(i) = 0.0;
				} else {
					edge_target_angle(i) = angles_json[i].get<double>() * M_PI / 180.0;
				}
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
				} else if (assignments[i] == "F" || assignments[i] == "U"){
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

		// fill in edge_adjacent_vertices using the help of faces_vertices and faces_edges
		if (params.contains("faces_vertices") && params.contains("faces_edges") && params.contains("edges_assignment")){
			std::vector<std::vector<double>> faces_vertices = params["faces_vertices"].template get<std::vector<std::vector<double>>>();
			std::vector<std::vector<double>> faces_edges = params["faces_edges"].template get<std::vector<std::vector<double>>>();
			std::vector<std::string> assignments = params["edges_assignment"].template get<std::vector<std::string>>();

			edge_adjacent_vertices.resize(E.rows(), 4);
			
			// Iterate through all edges and for every non-border-edge figure out it's 4 defining vertices
			for (int currentEdgeID = 0; currentEdgeID < E.rows(); currentEdgeID++){
				if (assignments[currentEdgeID] == "B"){
					edge_adjacent_vertices.row(currentEdgeID) << -1, -1, -1, -1;
				} else {
					int first_face_ID = -1;
					int second_face_ID = -1;

					// Iterate through all faces to find the (exactly) 2 faces which contain the current edge id
					for (int facesIterator = 0; facesIterator < faces_edges.size(); facesIterator++){
						if (std::find(faces_edges[facesIterator].begin(), faces_edges[facesIterator].end(), currentEdgeID) != faces_edges[facesIterator].end()){
							if (first_face_ID == -1){
								first_face_ID = facesIterator;
							} else {
								second_face_ID = facesIterator;
							}
						}
					}

					std::vector<double> first_face = faces_vertices[first_face_ID];
					std::vector<double> second_face = faces_vertices[second_face_ID];

					// Get start and end vertex of the edge
					int v_start = E(currentEdgeID, 0);
					int v_end = E(currentEdgeID, 1);

					// Helper to get third vertex not part of edge
					auto find_third_vertex = [&](const std::vector<double>& face) -> int {
						for (double v : face) {
							int vi = (int)v;
							if (vi != v_start && vi != v_end) return vi;
						}
						return -1; // Should not happen
					};

					int third1 = find_third_vertex(first_face);
					int third2 = find_third_vertex(second_face);

					Eigen::Vector3d a = V.row(v_start);
					Eigen::Vector3d b = V.row(v_end);
					Eigen::Vector3d c1 = V.row(third1);
					Eigen::Vector3d c2 = V.row(third2);

					// Calculate normals to triangle (v_start, v_end, v_third)
					Eigen::Vector3d edge_vec = b - a;
					Eigen::Vector3d n1 = edge_vec.cross(c1 - a);
					Eigen::Vector3d n2 = edge_vec.cross(c2 - a);

					// Check which of them points up / down
					Eigen::Vector3d Z_axis_normal(0, 0, 1);
					bool n1_pos = n1.dot(Z_axis_normal) > 0;
					bool n2_pos = n2.dot(Z_axis_normal) > 0;

					int v_left, v_right;

					// Assign left and right depending on up and down
					if (n1_pos && !n2_pos) {
						v_left = third1;
						v_right = third2;
					} else if (!n1_pos && n2_pos) {
						v_left = third2;
						v_right = third1;
					} else {
						std::cout << "something went wrong trying to deduce edce adjacent vertex ordering, defaulting..." << std::endl;
						v_left = third1;
						v_right = third2;
					}

					edge_adjacent_vertices.row(currentEdgeID) << v_right, v_left, v_start, v_end;	
				}
			}
		} else {
			std::cout << "Missing parameter \"faces_vertices\" or \"faces_edges\" or \"edges_assignment\"..." << std::endl;
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

		// Add the edge adjacent vertices, used in the calculation for the crease constraints
		edge_adjacent_vertices.resize(5, 4);
		edge_adjacent_vertices.row(0) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(1) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(2) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(3) << -1, -1, -1, -1;
		edge_adjacent_vertices.row(4) << 0, 3, 2, 1;
	}

}