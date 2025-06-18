#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <json.hpp>
#include <fstream>

// Structs to keep track of simulation state and parameters
struct SimulationData {
    double t = 0;                                   // Simulation time
    Eigen::VectorXd q;                              // Generalized Vertex Coordinates, size 3*x vector
    Eigen::VectorXd qdot;                           // Generalized Vertex Velocities
    Eigen::SparseMatrix<double> M;                  // Sparse Mass Matrix.
    Eigen::MatrixXd V;                              // Vertices of the CP, nx3 matrix
    Eigen::MatrixXi F;                              // Faces of the CP, mx3 matrix
    Eigen::MatrixXi E;                              // Edges of the CP (Springs), ex2 matrix
    Eigen::VectorXd edge_target_angle;              // Final Target fold angle for each edge
    Eigen::VectorXd curr_theta;                     // Current Target fold angle for each edge (CURRENTLY UNUSED)
    Eigen::VectorXd l0;                             // Original length of the Edges, size e vector
    Eigen::MatrixXd alpha0;                         // Nominal angles in flat state of CP, in same order as face vertices
    Eigen::VectorXd k_axial;                        // Per edge stiffness constant
    Eigen::VectorXd k_crease;                       // Per crease stiffness constant
    Eigen::MatrixXi edge_adjacent_vertices;         // For each edge, stores the four vertices making up the two triangles which meet at the edge. The order is: Right vertex, Left Vertex, Start Vertex, End Vertex. It is {-1, -1, -1, -1} for border edges
    Eigen::MatrixXi face_adjacent_edges;            // For each face stores the three face adjacent edges as provided by the .fold field faces_edges
};

struct SimulationParams {
    // This controls the physics part of the simulation
    double dt;                  // Time Step
    double vertexMass;          // Per vertex mass. Currently constant
    double EA;                  // Axial stiffness parameter, used in calculating axial stiffness
    double k_fold;              // Stiffness for a mountain or valley crease (Should be much smaller than the axial stiffness)
    double k_facet;             // Stiffness for a facet crease
    double k_face;              // Stiffness for the face constraints
    double zeta;                // Parameter in the damping ratio from the paper
    Eigen::Vector3d g;          // Gravity force vector

    // This controls broader simulation parameters. Self explanatory (I hope)
    bool ENABLE_STRAIN_VISUALIZATION;
    std::string STRAIN_TYPE;
    bool ENABLE_DYNAMIC_SIMULATION;
    bool ENABLE_GRAVITY;
    bool USE_IMPLICIT_EULER;
};

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

