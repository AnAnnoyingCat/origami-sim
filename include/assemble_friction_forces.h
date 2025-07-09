#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <axial_constraints.h>
#include <iostream>
#include <parameters.h>
#include <ipc/ipc.hpp>
#include <ipc/potentials/friction_potential.hpp>


/**
 * @brief Calculates size 3*n vector of forces acting from collisions with the ground using IPC-Toolkit. Assumes collisions were previously built
 * 
 */
void assemble_friction_forces_IPC(Eigen::VectorXd &f, SimulationParams& simulationParams, SimulationData& simulationData);