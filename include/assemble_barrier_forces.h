#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <axial_constraints.h>
#include <iostream>
#include <parameters.h>
#include <ipc/ipc.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/barrier/adaptive_stiffness.hpp>
#include <ipc/utils/world_bbox_diagonal_length.hpp>

/**
 * @brief Calculate the size 3*n vector of forces acting from collision with the ground
 * 
 * @param f 
 * @param q 
 * @param min_barrier_distance 
 */
void assemble_ground_barrier_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, double min_barrier_distance, SimulationParams& simulationParams, SimulationData simulationData);

void get_barrier_force_for_vertex(Eigen::Vector3d &f, const Eigen::Vector3d q0, double d, double k_barrier);

/**
 * @brief Calculates size 3*n vector of forces acting from collisions with the ground using IPC-Toolkit. Assumes collisions were previously built
 * 
 */
void assemble_barier_forces_IPC(Eigen::VectorXd &f, SimulationParams& simulationParams, SimulationData& simulationData, bool first_time);