#include <dynamic_target_angle.h>
#include <map>

// std::map<Interval, std::map<int, double>> time_angles{
//     {Interval(0.00, 0.17), {{1, 90}, {3, -80}}}
// };
std::vector<FoldInstruction> fold_timeline;
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

void setupTimeDependantTargetAngle(std::vector<FoldInstruction> in, double numEdges){
	fold_timeline = in;
	per_fold_position_in_timeline.resize(numEdges);
	for (int i = 0; i < numEdges; i++){
		// At the beginning none of the edges are in the timeline
		per_fold_position_in_timeline[i] = find_pos_in_timeline(i, 0);
	}
}

void calculateDynamicTargetAngle(Eigen::VectorXd &edge_target_angle, double t) {
	// Iterate through every edge, check which interval it's on and lerp the angle based on the interval and current time
	for (int currEdge = 0; currEdge < edge_target_angle.size(); currEdge++){
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
		}

		FoldInstruction currentInstructions = fold_timeline[per_fold_position_in_timeline[currEdge]];
		// Linearly interpolate currEdge target angle based on currentInstructions0

	}
}