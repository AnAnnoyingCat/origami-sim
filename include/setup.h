#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <json.hpp>
#include <fstream>
#include <parameters.h>
#include <igl/opengl/glfw/Viewer.h>

/**
 * @brief Read all the simulation parameters from file filename
 * 
 * @param filename          Path to the file to be read
 * @param simulationParams  Struct in which the params will be returned
 */
void setup_simulation_params(std::string filename, SimulationParams& simulationParams);

/**
 * @brief Read all the information for the CP at filename and load it into the simulation, precompute some values
 * 
 * @param filename 
 * @param simulationParams 
 * @param simulationData 
 */
void setup_mesh(std::string filename, SimulationParams& simulationParams, SimulationData& simulationData);

// NaN -> Free swinging crease

/**
 * @brief Read all the time dependent target angle information and set up a fold timeline to be used later
 * 
 * @param filename      
 * @param edge_target_angle 
 */
void setup_dynamic_target_angles(std::string filename, Eigen::VectorXd& edge_target_angle);

/**
 * @brief Sets up a large floor at z=0
 * 
 */
int setup_floor(igl::opengl::glfw::Viewer& viewer);
