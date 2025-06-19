#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>

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
    double gluefactor;          // Glued edge i will have stiffness edge_stiffness(i) * gluefactor
    double zeta;                // Parameter in the damping ratio from the paper
    double min_barrier_distance;// Parameter which determines how close origami can get to the gorund before barrier forces engage
    double spawn_height;        // Controls how high up the flat paper spawns before falling to the ground. 
    Eigen::Vector3d g;          // Gravity force vector

    // This controls broader simulation parameters. Self explanatory (I hope)
    bool ENABLE_STRAIN_VISUALIZATION;
    std::string STRAIN_TYPE;
    bool ENABLE_DYNAMIC_SIMULATION;
    bool ENABLE_GRAVITY;
    bool USE_IMPLICIT_EULER;
    bool LOG_SIMULATION_TIME;
    bool USE_SNAPPING_GLUE_MODE;
};