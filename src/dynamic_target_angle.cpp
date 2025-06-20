#include <dynamic_target_angle.h>

// Vector of FoldInstructions containing all the keyframes for the different creases. This is ordered in a per crease fasion
std::vector<FoldInstruction> fold_timeline;

// Each fold's current position in the timeline.
std::vector<int> per_fold_position_in_timeline;

// Finds the next instruction for edgeNr in fold_timeline starting at startPos, returns -1 if not found
int find_pos_in_timeline(int edgeNr, int startPos){
	for (int i = startPos; i < fold_timeline.size(); i++){
		if (fold_timeline[i].fold_number == edgeNr){
			return i;
		}
	}
	return -1;
}

void setupFoldTimeline(std::vector<FoldInstruction> in, double numEdges){
	fold_timeline = in;
	per_fold_position_in_timeline.resize(numEdges);
	for (int i = 0; i < numEdges; i++){
		// At the beginning none of the edges are in the timeline
		per_fold_position_in_timeline[i] = find_pos_in_timeline(i, 0);
	}
}

void calculateDynamicTargetAngle(SimulationData& simulationData, SimulationParams& simulationParams) {
	double t = simulationData.t;
	Eigen::VectorXd q = simulationData.q;
	Eigen::MatrixXi edge_adjacent_vertices = simulationData.edge_adjacent_vertices;
	// Iterate through every edge, check which interval it's on and lerp the angle based on the interval and current time
	for (int currEdge = 0; currEdge < simulationData.edge_target_angle.size(); currEdge++){
		// Check if (any more) instructions for current edge exist
		if (per_fold_position_in_timeline[currEdge] == -1){
			continue;
		}
		
		// Check if current edge left its interval and interval needs updating
		if (!fold_timeline[per_fold_position_in_timeline[currEdge]].time.contains(t)){
			per_fold_position_in_timeline[currEdge] = find_pos_in_timeline(currEdge, per_fold_position_in_timeline[currEdge] + 1);
			if (per_fold_position_in_timeline[currEdge] == -1){
				// That was the last instruction
				continue;
			}
			simulationData.k_axial(currEdge) = simulationParams.EA / simulationData.l0(currEdge);
			//std::cout << "Fold " << currEdge << " got new instructions. Reset stiffness to: " << simulationData.k_axial(currEdge) << std::endl;
		}
		
		FoldInstruction currentInstructions = fold_timeline[per_fold_position_in_timeline[currEdge]];
		Interval time = currentInstructions.time;      
		double start_angle = currentInstructions.start_angle; 
		double end_angle = currentInstructions.end_angle;

		if (currentInstructions.mode.has_value()){
			std::cout << "we have a mode: " << currentInstructions.mode.value() << std::endl;
			if (currentInstructions.mode == "glued"){
				// Fold is now glued. 

				if (simulationParams.USE_SNAPPING_GLUE_MODE){
					// Set the target angles to be the percise expected target angle (what we *should* have reached up to here)
					fold_timeline[per_fold_position_in_timeline[currEdge]].end_angle = start_angle;
					end_angle = start_angle;
				} else {
					// Calculate what angle we actually currently have and glue that angle to be fix
					Eigen::Vector3d n1, n2;
					Eigen::Vector3d q1 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 0));
					Eigen::Vector3d q2 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 1));
					Eigen::Vector3d q3 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 2));
					Eigen::Vector3d q4 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 3));
					getNormal(n1, q1, q4, q3);
					getNormal(n2, q2, q3, q4);
					Eigen::Vector3d crease_dir = (q4 - q3).normalized();
					double current_theta = std::atan2((n1.cross(n2)).dot(crease_dir), n1.dot(n2));

					fold_timeline[per_fold_position_in_timeline[currEdge]].start_angle = current_theta * 180 / M_PI;
					fold_timeline[per_fold_position_in_timeline[currEdge]].end_angle = current_theta * 180 / M_PI;
					start_angle = fold_timeline[per_fold_position_in_timeline[currEdge]].start_angle;
					end_angle = start_angle;
				}
				
				// Update the fold stiffness to the glued stiffness
				simulationData.k_crease(currEdge) = simulationParams.k_fold * simulationParams.gluefactor;

				// Remove the glued keyword so this code only gets applied once (which it needs to)
				fold_timeline[per_fold_position_in_timeline[currEdge]].mode.reset();
				
			} else if (currentInstructions.mode == "free"){
				// Freely swinging creases have a target angle of NaN, which nicely propagates through all equations
				fold_timeline[per_fold_position_in_timeline[currEdge]].end_angle = nan("");
			} else if (currentInstructions.mode == "toggle_gravity") {
				// Toggle gravity and reset instructions
				std::cout << "Toggling gravity!" << std::endl;
				simulationParams.ENABLE_GRAVITY = !simulationParams.ENABLE_GRAVITY;
				fold_timeline[per_fold_position_in_timeline[currEdge]].mode.reset();
			} else {
				std::cerr << "unknown fold mode in fold instruction" << std::endl;
			}
		}

		if (std::isnan(start_angle) && !std::isnan(end_angle)){
			// Fold stops swinging freely and should lerp to target angle
			// Calculate current edge angle -> this will be the start angle
			Eigen::Vector3d n1, n2;
			Eigen::Vector3d q1 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 0));
			Eigen::Vector3d q2 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 1));
			Eigen::Vector3d q3 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 2));
			Eigen::Vector3d q4 = q.segment<3>(3 * edge_adjacent_vertices(currEdge, 3));
			getNormal(n1, q1, q4, q3);
			getNormal(n2, q2, q3, q4);
			Eigen::Vector3d crease_dir = (q4 - q3).normalized();
			double current_theta = std::atan2((n1.cross(n2)).dot(crease_dir), n1.dot(n2));

			fold_timeline[per_fold_position_in_timeline[currEdge]].start_angle = current_theta * 180 / M_PI;
			start_angle = current_theta * 180 / M_PI;

		} else if (std::isnan(end_angle) && !std::isnan(start_angle)){
			// Fold starts to swing freely and should freely swing for the entire interval
			// This should happen automatically since NaN will propagate through the following calculation and thus edge_target_angle(currEdge) should end up being NaN
			//std::cout << "This fold is now a NaN: " << currEdge << std::endl;
		}

		// lerp the angle based on the time
		double target_in_degrees = start_angle + (end_angle - start_angle) * (std::min((t - time.begin) / ((time.end - time.begin)* 0.8), 1.0));

		simulationData.edge_target_angle(currEdge) = target_in_degrees * M_PI / 180;
	}
}