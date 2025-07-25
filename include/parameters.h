#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <ipc/ipc.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>

// Structs to keep track of simulation state and parameters
struct SimulationData {

    // Model Variables
    Eigen::VectorXd q;                              // Generalized Vertex Coordinates, size 3*x vector
    Eigen::VectorXd qdot;                           // Generalized Vertex Velocities
    Eigen::SparseMatrix<double> M;                  // Sparse Mass Matrix.
    Eigen::MatrixXd V;                              // Vertices of the CP, nx3 matrix
    Eigen::MatrixXi F;                              // Faces of the CP, mx3 matrix
    Eigen::MatrixXi E;                              // Edges of the CP (Springs), ex2 matrix

    // Ground Variables
    Eigen::MatrixXd ground_V;                       // Ground vertex positions
    Eigen::MatrixXi ground_F;                       // Ground faces
    Eigen::MatrixXi ground_E;                       // Ground Edges

    // IPC Variables
    Eigen::MatrixXd rest_positions;                 // Undeformed vertex positions for IPC <- Includes both model and ground. always in the order << V, ground_V
    Eigen::MatrixXd deformed_vertices;              // Deformed vertex positions for IPC <- both model and ground
    ipc::CollisionMesh collision_mesh;              // Collision mesh used by IPC
    ipc::Collisions collisions;                     // IPCs active collisions
    ipc::BarrierPotential barrier_potential;        // IPC's barrier potential used in the whole simulation  (Initialized later)
    ipc::FrictionPotential friction_potential;      // IPC's friction potential used in the whole simulation (Initialized later)
    ipc::FrictionCollisions friction_collisions;    // Active collisions used by IPC for calculating friction
    double bbox_diagonal;                           // Lets IPC track the scale of the simulation. call this with only the simulation mesh because the floor is overly large
    double max_barrier_stiffness;                   // Initial guess of barrier stiffness and max it can reach this newton solve
    double barrier_stiffness;                       // Current barrierr stiffness used by my solver this iteration
    double prev_distance;                           // Previous minimal distance between any primitives
    double curr_distance;                           // Current minimal distance between any primitives

    
    // Simulation Variables
    double t = 0;                                   // Simulation time
    Eigen::VectorXd edge_target_angle;              // Final Target fold angle for each edge
    Eigen::VectorXd curr_theta;                     // Current Target fold angle for each edge (CURRENTLY UNUSED)
    Eigen::VectorXd l0;                             // Original length of the Edges, size e vector
    Eigen::MatrixXd alpha0;                         // Nominal angles in flat state of CP, in same order as face vertices
    Eigen::VectorXd k_axial;                        // Per edge stiffness constant
    Eigen::VectorXd k_crease;                       // Per crease stiffness constant
    Eigen::MatrixXi edge_adjacent_vertices;         // For each edge, stores the four vertices making up the two triangles which meet at the edge. The order is: Right vertex, Left Vertex, Start Vertex, End Vertex. It is {-1, -1, -1, -1} for border edges
    Eigen::MatrixXi face_adjacent_edges;            // For each face stores the three face adjacent edges as provided by the .fold field faces_edges

    // Constructor
    SimulationData()
        : barrier_potential(1.0),
          friction_potential(1.0)
    {}
};

struct SimulationParams {
    // This controls the physics part of the simulation
    double dt;                  // Time Step
    double vertexMass;          // Per vertex mass. Currently constant
    double EA;                  // Axial stiffness parameter, used in calculating axial stiffness
    double k_fold;              // Stiffness for a mountain or valley crease (Should be much smaller than the axial stiffness)
    double k_facet;             // Stiffness for a facet crease
    double k_face;              // Stiffness for the face constraints
    double k_barrier;           // Stiffness for the barrier force
    double mu;                  // Global coefficient of friction
    double eps_v;               // Friction threshold in units of velocity
    double gluefactor;          // Glued edge i will have stiffness edge_stiffness(i) * gluefactor
    double zeta;                // Parameter in the damping ratio from the paper
    double min_barrier_distance;// Parameter which determines how close origami can get to the gorund before barrier forces engage
    double spawn_height;        // Controls how high up the flat paper spawns before falling to the ground. 
    Eigen::Vector3d g;          // Gravity force vector
    bool simulating;            // If this is ever set to false, all threads instantly terminates
    double sim_zoom_level;      // The starting zoom level of the simulation
    bool enable_floor;          // Toggles whether or not the floor is visible
    bool enable_barrier;        // Toggles whether or not to use barrier forces
    bool enable_friction;       // Toggles whether or not to use friction forces
    bool enable_auto_k_barrier; // Toggles whether or not to use IPC's auto barrier stiffness
    bool loop_timeline;         // Toggles whether or not the timeline will loop, jumping to t=0 after reaching last instruction
    bool center_mesh;           // Toggles whether or not the model should be centered every frame to prevent drifting

    // This controls broader simulation parameters. Self explanatory (I hope)
    bool ENABLE_STRAIN_VISUALIZATION;
    std::string STRAIN_TYPE;
    bool ENABLE_DYNAMIC_SIMULATION;
    bool ENABLE_GRAVITY;
    bool USE_IMPLICIT_EULER;
    bool LOG_SIMULATION_TIME;
    bool USE_SNAPPING_GLUE_MODE;
    bool LOG_FORCES;

};