#include <setup.h>
#include <iostream>
#include <trig_helper_functions.h>
#include <dynamic_target_angle.h>
#include <vector>

using json = nlohmann::json;

void setup_simulation_params(std::string filename, SimulationParams& simulationParams){
	// Use nlohmann JSON to grab the simulation parameters from the specified filename.
	std::ifstream file(filename);
	if (file){
		json params = json::parse(file);

		if (params.contains("dt")){
			simulationParams.dt = params["dt"].template get<double>();
		} else {
			simulationParams.dt = 0.001;
		}
		if (params.contains("vertexMass")){
			simulationParams.vertexMass = params["vertexMass"].template get<double>();
		} else {
			simulationParams.vertexMass = 1;
		}
		if (params.contains("EA")){
			simulationParams.EA = params["EA"].template get<double>();
		} else {
			simulationParams.EA = 5e4;
		}
		if (params.contains("k_fold")){
			simulationParams.k_fold = params["k_fold"].template get<double>();
		} else {
			simulationParams.k_fold = 1e3;
		}
		if (params.contains("k_facet")){
			simulationParams.k_facet = params["k_facet"].template get<double>();
		} else {
			simulationParams.k_facet = 1e3;
		}
		if (params.contains("k_face")){
			simulationParams.k_face = params["k_face"].template get<double>();
		} else {
			simulationParams.k_face = 2e2;
		}
		if (params.contains("zeta")){
			simulationParams.zeta = params["zeta"].template get<double>();
		} else {
			simulationParams.zeta = 0.25;
		}
		if (params.contains("min_barrier_distance")){
			simulationParams.min_barrier_distance = params["min_barrier_distance"].template get<double>();
		} else {
			simulationParams.min_barrier_distance = 0.1;
		}
		if (params.contains("visualize_strain")){
			simulationParams.ENABLE_STRAIN_VISUALIZATION = params["visualize_strain"].template get<bool>();
		} else {
			simulationParams.ENABLE_STRAIN_VISUALIZATION = false;
		}
		if (params.contains("strain_type")){
			simulationParams.STRAIN_TYPE = params["strain_type"].template get<std::string>();
		} else {
			simulationParams.STRAIN_TYPE = "face";
		}
		if (params.contains("enable_dynamic_simulation")){
			simulationParams.ENABLE_DYNAMIC_SIMULATION = params["enable_dynamic_simulation"].template get<bool>();
		} else {
			simulationParams.ENABLE_DYNAMIC_SIMULATION = true;
		}
		if (params.contains("enable_gravity")){
			simulationParams.ENABLE_GRAVITY = params["enable_gravity"].template get<bool>();
		} else {
			simulationParams.ENABLE_GRAVITY = true;
		}
		if (params.contains("gravity")){
			for (int i = 0; i < 3; ++i){
				simulationParams.g(i) = params["gravity"][i].get<double>();
			}
		} else {
			simulationParams.g = Eigen::Vector3d(0, 0, -9.81); // default gravity
		}
		if (params.contains("use_implicit_euler")){
			simulationParams.USE_IMPLICIT_EULER = params["use_implicit_euler"].template get<bool>();
		} else {
			simulationParams.USE_IMPLICIT_EULER = true;
		}
		if (params.contains("glue_factor")){
			simulationParams.gluefactor = params["glue_factor"].template get<double>();
		} else {
			simulationParams.gluefactor = 1;
		}
		if (params.contains("log_simulation_time")){
			simulationParams.LOG_SIMULATION_TIME = params["log_simulation_time"].template get<bool>();
		} else {
			simulationParams.LOG_SIMULATION_TIME = false;
		}
		if (params.contains("use_snapping_glue_mode")){
			simulationParams.USE_SNAPPING_GLUE_MODE = params["use_snapping_glue_mode"].template get<bool>();
		} else {
			simulationParams.USE_SNAPPING_GLUE_MODE = true;
		}
		if (params.contains("spawn_height")){
			simulationParams.spawn_height = params["spawn_height"].template get<double>();
		} else {
			simulationParams.spawn_height = 1;
		}
		if (params.contains("log_forces")){
			simulationParams.LOG_FORCES = params["log_forces"].template get<bool>();
		} else {
			simulationParams.LOG_FORCES = false;
		}
		if (params.contains("show_floor")){
			simulationParams.show_floor = params["show_floor"].template get<bool>();
		} else {
			simulationParams.show_floor = false;
		}
		if (params.contains("k_barrier")){
			simulationParams.k_barrier = params["k_barrier"].template get<double>();
		} else {
			simulationParams.k_barrier = 1;
		}
		if (params.contains("sim_zoom_level")){
			simulationParams.sim_zoom_level = params["sim_zoom_level"].template get<double>();
		} else {
			simulationParams.sim_zoom_level = 250;
		}


		simulationParams.simulating = true;
		
	} else {
		std::cout << "error: params file not found" << std::endl;     
	}
}

void setup_mesh(std::string filename, SimulationParams& simulationParams, SimulationData& simulationData){
	std::ifstream file(filename);

	if (file){
		json params = json::parse(file);
		// Check for consistency of edge fields 
		if (params.contains("edges_vertices") && params.contains("edges_assignment") && params.contains("edges_foldAngle")){
			const auto& edges_vertices = params["edges_vertices"];
			const auto& edges_assignment = params["edges_assignment"];
			const auto& edges_foldAngle = params["edges_foldAngle"];
			
			if (!(edges_vertices.size() == edges_assignment.size() &&
			edges_vertices.size() == edges_foldAngle.size())) {
			std::cerr << "[ERROR] Mismatched edge data: "
					<< "vertices=" << edges_vertices.size()
					<< ", assignment=" << edges_assignment.size()
					<< ", angle=" << edges_foldAngle.size()
					<< std::endl;
			exit(EXIT_FAILURE);
			}
		}

		// Assign Vertices
		if (params.contains("vertices_coords")){
			std::vector<std::vector<double>> coords = params["vertices_coords"].template get<std::vector<std::vector<double>>>();

			simulationData.V.resize(coords.size(), 3);

			double largestCoord;
			// find largest coordinate and scale entire model to make largets coordinate to be around the ballpark of 100
			for (int i = 0; i < coords.size(); i++){
				std::cout << "simulationData.V(i, 0): " << coords[i][0] << "largestcoord: " << largestCoord << std::endl;
				largestCoord = largestCoord < coords[i][0] ? coords[i][0] : largestCoord;
				std::cout << "simulationData.V(i, 1): " << coords[i][1] << "largestcoord: " << largestCoord << std::endl;
				largestCoord = largestCoord < coords[i][1] ? coords[i][1] : largestCoord;
			}
			double scalefactor = 100 / largestCoord;
			std::cout << "largest coord: " << largestCoord << ", scale factor: " << scalefactor << std::endl;

			for (int i = 0; i < coords.size(); i++){
				simulationData.V(i, 0) = coords[i][0] * scalefactor;
				simulationData.V(i, 1) = coords[i][1] * scalefactor;
				if (coords[i].size() == 3){
					simulationData.V(i, 2) = coords[i][2] * scalefactor;
				} else {
					simulationData.V(i, 2) = simulationParams.spawn_height;
				}
			}
		} else {
			std::cout << "Missing parameter \"vertices_coords\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}

		// Assign faces
		if (params.contains("faces_vertices")){
			std::vector<std::vector<double>> faces = params["faces_vertices"].template get<std::vector<std::vector<double>>>();
			simulationData.F.resize(faces.size(), 3);
			for (int i = 0; i < faces.size(); i++){
				simulationData.F(i, 0) = faces[i][0];
				simulationData.F(i, 1) = faces[i][1];
				simulationData.F(i, 2) = faces[i][2];
			}
		} else {
			std::cout << "Missing parameter \"faces_vertices\"" << std::endl;
			std::cout << "we do have vertices coords. proof: " << simulationData.V << std::endl;
			throw std::runtime_error("Oh no...");
		}
		
		// Assign resting face angles
		simulationData.alpha0.resize(simulationData.F.rows(), 3);
		for (int i = 0; i < simulationData.F.rows(); i++){
			// Angle 1 of ith face is the angle at the first vertex, between the second and third
			getAngle(simulationData.alpha0(i, 0), simulationData.V.row(simulationData.F(i, 0)), simulationData.V.row(simulationData.F(i, 1)), simulationData.V.row(simulationData.F(i, 2)));
			// Angle 2 is at the second vertex, between the third and the first
			getAngle(simulationData.alpha0(i, 1), simulationData.V.row(simulationData.F(i, 1)), simulationData.V.row(simulationData.F(i, 2)), simulationData.V.row(simulationData.F(i, 0)));
			// Angle 3 is at the third vertex between the first and the second
			getAngle(simulationData.alpha0(i, 2), simulationData.V.row(simulationData.F(i, 2)), simulationData.V.row(simulationData.F(i, 0)), simulationData.V.row(simulationData.F(i, 1)));
		}

		// Assign edges
		if (params.contains("edges_vertices")){
			std::vector<std::vector<double>> edges = params["edges_vertices"].template get<std::vector<std::vector<double>>>();
			simulationData.E.resize(edges.size(), 2);
			for (int i = 0; i < edges.size(); i++){
				simulationData.E(i, 0) = edges[i][0];
				simulationData.E(i, 1) = edges[i][1];
			}
		} else {
			std::cout << "Missing parameter \"edges_vertices\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}

		// Set up faces → edges
		if (params.contains("faces_edges")){
			std::vector<std::vector<double>> faces_edges = params["faces_edges"].template get<std::vector<std::vector<double>>>();
			simulationData.face_adjacent_edges.resize(faces_edges.size(), 3);
			for (int i = 0; i < faces_edges.size(); i++){
				simulationData.face_adjacent_edges(i, 0) = faces_edges[i][0];
				simulationData.face_adjacent_edges(i, 1) = faces_edges[i][1];
				simulationData.face_adjacent_edges(i, 2) = faces_edges[i][2];
			}
		} else {
			std::cout << "Missing parameter \"edges_vertices\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}


		// Assign edge target angles. If edge doesn't have a target angle the angle is set to nan.
		if (params.contains("edges_foldAngle") && params.contains("edges_assignment")){
			const auto& angles_json = params["edges_foldAngle"];
			const auto& assignments = params["edges_assignment"];

			simulationData.edge_target_angle.resize(angles_json.size());

			for (int i = 0; i < angles_json.size(); i++){
				if (angles_json[i].is_null()){
					// null angle usually means border edge or facet crease. in any case, set 0.
					simulationData.edge_target_angle(i) = nan("");
				} else {
					if (assignments[i].get<std::string>() == "U" || assignments[i].get<std::string>() == "B"){
						simulationData.edge_target_angle(i) = nan("");
					} else {
						simulationData.edge_target_angle(i) = angles_json[i].get<double>() * M_PI / 180.0;
					}
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
			simulationData.l0.resize(simulationData.E.rows());
			simulationData.k_axial.resize(simulationData.E.rows());
			simulationData.k_crease.resize(simulationData.E.rows());
			for (int i = 0; i < simulationData.E.rows(); i++){
				int v0 = simulationData.E(i, 0);
				int v1 = simulationData.E(i, 1);
				Eigen::Vector3d difference = simulationData.V.row(v0) - simulationData.V.row(v1);

				// Set up l0
				simulationData.l0(i) = difference.norm();

				// Set up axial stiffness dependent on l0
				simulationData.k_axial(i) = simulationParams.EA / simulationData.l0(i);

				// Set up crease stiffness dependent on l0
				if (assignments[i] == "B"){
					simulationData.k_crease(i) = -1;
				} else if (assignments[i] == "F"){
					simulationData.k_crease(i) = simulationParams.k_facet * simulationData.l0(i);
				} else { // Edge type == "M" or "V"
					simulationData.k_crease(i) = simulationParams.k_fold * simulationData.l0(i);
				}
			}
			
		} else {
			std::cout << "Missing parameter \"edges_assignment\"" << std::endl;
			throw std::runtime_error("Oh no...");
		}
		

		// Initialize q and qdot
		simulationData.q.resize(simulationData.V.rows() * simulationData.V.cols());
		simulationData.qdot.resize(simulationData.V.rows() * simulationData.V.cols());

		Eigen::MatrixXd Vt = simulationData.V.transpose();
		simulationData.q = Eigen::Map<Eigen::VectorXd>(Vt.data(), Vt.rows() * Vt.cols());
		simulationData.qdot.setZero();

		// fill in edge_adjacent_vertices using the help of faces_vertices and faces_edges
		if (params.contains("faces_vertices") && params.contains("faces_edges") && params.contains("edges_assignment")){
			std::vector<std::vector<double>> faces_vertices = params["faces_vertices"].template get<std::vector<std::vector<double>>>();
			std::vector<std::vector<double>> faces_edges = params["faces_edges"].template get<std::vector<std::vector<double>>>();
			std::vector<std::string> assignments = params["edges_assignment"].template get<std::vector<std::string>>();

			simulationData.edge_adjacent_vertices.resize(simulationData.E.rows(), 4);
			
			// Iterate through all edges and for every non-border-edge figure out it's 4 defining vertices
			for (int currentEdgeID = 0; currentEdgeID < simulationData.E.rows(); currentEdgeID++){
				if (assignments[currentEdgeID] == "B"){
					simulationData.edge_adjacent_vertices.row(currentEdgeID) << -1, -1, -1, -1;
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
					int v_start = simulationData.E(currentEdgeID, 0);
					int v_end = simulationData.E(currentEdgeID, 1);

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

					Eigen::Vector3d a = simulationData.V.row(v_start);
					Eigen::Vector3d b = simulationData.V.row(v_end);
					Eigen::Vector3d c1 = simulationData.V.row(third1);
					Eigen::Vector3d c2 = simulationData.V.row(third2);

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

					simulationData.edge_adjacent_vertices.row(currentEdgeID) << v_right, v_left, v_start, v_end;	
				}
			}
		} else {
			std::cout << "Missing parameter \"faces_vertices\" or \"faces_edges\" or \"edges_assignment\"..." << std::endl;
			throw std::runtime_error("Oh no...");
		}



	} else {
		std::cerr << "There was an error reading your .fold file!" << std::endl;
	}
}

void setup_dynamic_target_angles(std::string filename, Eigen::VectorXd& edge_target_angle){
	std::ifstream file(filename);
	// save all the fold information in a custom struct
	std::vector<FoldInstruction> foldTimeline;

	if (file){
		json params = json::parse(file);
		//std::cout << "Found file! " << std::endl;
		std::map<int, double> last_timestamp;

		// Store the last angle of each crease for future use
		std::map<int, double> lastAngle;
		for (int i = 0; i < edge_target_angle.size(); i++){
			lastAngle[i] = 0.0;
			last_timestamp[i] = 0.0;
		}

		// Store the keys as doubles so they're sorted correctly
		std::map<double, nlohmann::json> sorted_params;

		for (auto& [key, value] : params.items()) {
			sorted_params[std::stod(key)] = value;
		}

		// Iterate through all the keyframes and add them to the datastructure
		for (auto& [timestamp, value] : sorted_params){
			//std::cout << "Key (double) : " << timestamp << std::endl;
			for (const auto& pair : value){
				int creaseNr = pair[0];

				double targetAngle;
				std::optional<std::string> mode;
				if (pair[1].is_null()){
					targetAngle = nan("");
				} else if (pair[1].is_number()) {
					targetAngle = pair[1].get<double>();
				} else if (pair[1].is_string()){
					std::string val = pair[1].get<std::string>();
					if (val == "glued" || val == "free" || val == "toggle_gravity") {
						targetAngle = nan("");
						mode = val;
					} else {
						std::cerr << "Warning: Unknown mode string '" << val << "' for fold " << creaseNr << "\n";
					}
				} else {
					std::cerr << "Warning: unknown fold instruction" << std::endl;
				}

				FoldInstruction instr = {
					Interval(last_timestamp[creaseNr], timestamp),
					creaseNr,
					lastAngle[creaseNr],
					targetAngle,
					mode
				};
				//std::cout << "  Pair: (" << creaseNr << ", " << targetAngle << ")" << std::endl;
				foldTimeline.push_back(instr);
				
				lastAngle[creaseNr] = targetAngle;
				last_timestamp[creaseNr] = timestamp;
			}
		}
	} else {
		std::cout << "No activation profile provided. Using default: Fully folded at t = 1" << std::endl;
		
		for (int i = 0; i < edge_target_angle.size(); i++){
			foldTimeline.push_back({Interval(0, 1), i, 0, edge_target_angle(i)});
		}
	}

	for (const auto& instruction : foldTimeline) {
        std::cout << "Fold " << instruction.fold_number 
                  << ": From " << instruction.start_angle << "° to " << instruction.end_angle << "°"
                  << " during [" << instruction.time.begin << " - " << instruction.time.end << "]\n";
    }

	setupFoldTimeline(foldTimeline, edge_target_angle.size());
}

int setup_floor(igl::opengl::glfw::Viewer& viewer) {
    // Create vertices for a square floor at z = 0
    Eigen::MatrixXd floor_vertices(4, 3);
    floor_vertices << -1000, -1000, 0,
                       1000, -1000, 0,
                       1000,  1000, 0,
                      -1000,  1000, 0;

    Eigen::MatrixXi floor_faces(2, 3);
    floor_faces << 0, 1, 2,
                   0, 2, 3;

    int floor_id = viewer.append_mesh();
    viewer.data(floor_id).set_mesh(floor_vertices, floor_faces);

	// Floor visual settings
    viewer.data(floor_id).set_colors(Eigen::RowVector3d(0.6, 0.6, 0.6));  // light grey
    viewer.data(floor_id).show_lines = false;
    viewer.data(floor_id).show_texture = false;
    viewer.data(floor_id).show_faces = true;
	return floor_id;
}