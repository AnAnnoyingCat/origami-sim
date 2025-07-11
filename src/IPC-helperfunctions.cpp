#include <IPC-helperfunctions.h>

void make_collision_mesh(SimulationData& simulationData, SimulationParams& simulationParams) {
	if (simulationParams.enable_floor){
		// Initialize the rest positions
		Eigen::MatrixXd rest_positions(simulationData.V.rows() + simulationData.ground_V.rows(), 3);
		rest_positions.topRows(simulationData.V.rows()) = simulationData.V;
		rest_positions.bottomRows(simulationData.ground_V.rows()) = simulationData.ground_V;
		simulationData.rest_positions = rest_positions;

		// Initialize total edges and faces
		Eigen::MatrixXi total_edges(simulationData.E.rows() + simulationData.ground_E.rows(), 2);
		total_edges.topRows(simulationData.E.rows()) = simulationData.E;
		Eigen::MatrixXi offset_ground_E = simulationData.ground_E.array() + simulationData.V.rows();
		total_edges.bottomRows(simulationData.ground_E.rows()) = offset_ground_E;

		Eigen::MatrixXi total_faces(simulationData.F.rows() + simulationData.ground_F.rows(), 3);
		total_faces.topRows(simulationData.F.rows()) = simulationData.F;
		Eigen::MatrixXi offset_ground_F = simulationData.ground_F.array() + simulationData.V.rows();
		total_faces.bottomRows(simulationData.ground_F.rows()) = offset_ground_F;
		
		ipc::CollisionMesh new_collision_mesh(rest_positions, total_edges, total_faces);
		simulationData.collision_mesh = new_collision_mesh;
	} else {
		simulationData.rest_positions = simulationData.V;
		ipc::CollisionMesh new_collision_mesh(simulationData.V, simulationData.E, simulationData.F);
		simulationData.collision_mesh = new_collision_mesh;
	}
}

void get_deformed_positions(SimulationData& simulationData, SimulationParams& simulationParams) {
	if (simulationParams.enable_floor){
		Eigen::MatrixXd deformed_vertices(simulationData.V.rows() + simulationData.ground_V.rows(), 3);
		deformed_vertices.topRows(simulationData.V.rows()) = simulationData.V;
		deformed_vertices.bottomRows(simulationData.ground_V.rows()) = simulationData.ground_V;
		simulationData.deformed_vertices = deformed_vertices;
	} else {
		simulationData.deformed_vertices = simulationData.V;
	}
}