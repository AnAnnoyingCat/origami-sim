#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <parameters.h>

/**
 * @brief Creates a collision mesh from the current Origami model and the floor
 * 
 */
void make_collision_mesh(SimulationData& simulationData, SimulationParams& simulationParams);

/**
 * @brief Returns the deformed Origami mesh with the floor appended to it
 * 
 * @param vertices 
 * @param simulationData 
 */
void get_deformed_positions(SimulationData& simulationData, SimulationParams& simulationParams);