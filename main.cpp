// C++ and libraries
#include <igl/opengl/glfw/Viewer.h>
#include <chrono>
#include <thread>
#include <Eigen/Dense>

// Project Functions
#include <setup.h>
#include <assemble_edge_forces.h>
#include <assemble_crease_forces.h>
#include <d2V_axial_dq2.h>
#include <axial_constraints.h>
#include <forward_Euler.h>
#include <fully_implicit_euler.h>
#include <assemble_edge_stiffness.h>
#include <make_mass_matrix.h>
#include <assemble_damping_forces.h>
#include <assemble_face_forces.h>
#include <trig_helper_functions.h>
#include <strain_calculations.h>
#include <dynamic_target_angle.h>
#include <assemble_crease_stiffness.h>
#include <assemble_damping_stiffness.h>
#include <assemble_gravity_forces.h>
#include <assemble_gravity_stiffness.h>

#include <finite_difference_tester.h>

// Simulation state
bool simulating = true;
Eigen::VectorXd q;                  // Generalized Vertex Coordinates, size 3*x vector
Eigen::VectorXd qdot;               // Generalized Vertex Velocities
Eigen::SparseMatrix<double> M;      // Sparse Mass Matrix.

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
Eigen::Vector3d g;                              // Gravity force vector

double t = 0;               // Simulation Time

double dt;                  // Time Step
double vertexMass;          // Per vertex mass. Currently constant
double EA;                  // Axial stiffness parameter, used in calculating axial stiffness
double k_fold;              // Stiffness for a mountain or valley crease (Should be much smaller than the axial stiffness)
double k_facet;             // Stiffness for a facet crease
double k_face;              // Stiffness for the face constraints
double zeta;                // Parameter in the damping ratio from the paper


// Working memory for integrator
Eigen::SparseMatrix<double> tmp_stiffness;
Eigen::VectorXd tmp_force;

// Pointer to the viewer for update reasons
igl::opengl::glfw::Viewer* viewer_ptr = nullptr;

// Debug flags
bool ENABLE_STRAIN_VISUALIZATION;
std::string STRAIN_TYPE;
bool ENABLE_DYNAMIC_SIMULATION;
bool ENABLE_GRAVITY;
bool USE_IMPLICIT_EULER;

/// @brief If enabled in config, visualizes strain in the mesh using l0 or alpha0
void visualizeStrain(){
    while(simulating){
        if (viewer_ptr){
            Eigen::MatrixXd C;
            if (STRAIN_TYPE == "face"){
                calculateFaceAngleStrain(C, F, q, alpha0);
            } else if (STRAIN_TYPE == "edge"){
                calculateAxialDeformationStrain(C, F, E, q, l0, face_adjacent_edges);
            }
            viewer_ptr->data().set_colors(C);
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
}

/// @brief Helper function which updates V using q
void updateV(Eigen::MatrixXd& V, Eigen::VectorXd q){
	int n = V.rows();  
	for (int i = 0; i < n; i++){
		V.row(i) = q.segment<3>(3 * i);
	}
}

/// @brief Move the center of mass back to the middle of the screen
void centerMesh(){ 
    Eigen::Vector3d center;
    center.setZero();
    for (int currVertex = 0; currVertex < V.rows(); currVertex++){
        center += q.segment<3>(3 * currVertex);
    }
    center /= V.rows();
    for (int currVertex = 0; currVertex < V.rows(); currVertex++){
        q.segment<3>(3 * currVertex) -= center;
    }
}

/// @brief Main simulation loop, running in a seperate thread
void simulate(){
    while (simulating){
        
        // If simulation type is dynamic, recalculate the target angle for the curren frame
        if (ENABLE_DYNAMIC_SIMULATION){
            calculateDynamicTargetAngle(edge_target_angle, t, q, edge_adjacent_vertices);
        }

        auto forces = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            // Set f to zero and then add all the forces to it
            f.resize(q.size());
            f.setZero();

            // Get the basic forces
            assemble_edge_forces(f, q, E, l0, k_axial);                                         // Keep edge lengths constant
            assemble_damping_forces(f, qdot, E, k_axial, zeta);                                 // Apply viscous dampening
            assemble_crease_forces(f, q, edge_adjacent_vertices, k_crease, edge_target_angle);  // Calculate crease driving force

            // If graivty is enabled, get that too
            if (ENABLE_GRAVITY){
                assemble_gravity_forces(f, g, vertexMass);
            }
        };

        auto stiffness = [&](Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            K.resize(q.size(), q.size());
            K.setZero();

            // Get all the basic stiffnessess
            assemble_edge_stiffness(K, q, V, E, l0, k_axial);
            assemble_crease_stiffness(K, q, edge_adjacent_vertices, k_crease, edge_target_angle);
            assemble_damping_stiffness(K, qdot, E, k_axial, zeta);

            /* If gravity is enabled, and uses a more complicated than normal computation, get that too
            if (ENABLE_GRAVITY){
                assemble_gravity_stiffness(K);
            }
            */
            
        };

        // Time integration
        if (USE_IMPLICIT_EULER){
            implicit_euler(q, qdot, dt, M, forces, stiffness, tmp_force, tmp_stiffness);
        } else {
            forward_euler(q, qdot, dt, forces, tmp_force);
        }
        
        

        // Make sure the floating mesh doesn't drift off with gravity disabled
        if (!ENABLE_GRAVITY){
            centerMesh();
        }

        // update vertex positions from q and increment time
        updateV(V, q);
        
        // Update the viewer's mesh
        if (viewer_ptr){
            viewer_ptr->data().set_vertices(V);
            viewer_ptr->data().compute_normals();
        }

        // Next time step
        t += dt;
        std::cout << "t: " << t << std::endl;
    }
}

int main(int argc, char *argv[])
{    

    // Read args into a vector
    std::vector<std::string> args;
    std::copy(argv + 1, argv + argc, std::back_inserter(args));

    // Set up parameters and read CP
    if (args.size() == 2){
        // Both CP and params provided
        setup_simulation_params(args[1], dt, vertexMass, EA, k_fold, k_facet, k_face, zeta, ENABLE_STRAIN_VISUALIZATION, STRAIN_TYPE, ENABLE_DYNAMIC_SIMULATION, ENABLE_GRAVITY, g, USE_IMPLICIT_EULER);
        setup_mesh(args[0], q, qdot, V, F, alpha0, E, edge_target_angle, l0, edge_adjacent_vertices, k_axial, k_crease, EA, k_fold, k_facet, k_face, face_adjacent_edges);

        // Get path to activation profiles
        args[0].replace(args[0].find("crease_patterns"), std::string("crease_patterns").length(), "activation_profiles");
        args[0].replace(args[0].find(".fold"), 5, ".json");
        setup_dynamic_target_angles(args[0], edge_target_angle);

    } else if (args.size() == 1){
        // Only crease pattern provided
        std::cout << "Using default parameters" << std::endl;
        setup_simulation_params("../data/simulation_params/default-params.json", dt, vertexMass, EA, k_fold, k_facet, k_face, zeta, ENABLE_STRAIN_VISUALIZATION, STRAIN_TYPE, ENABLE_DYNAMIC_SIMULATION, ENABLE_GRAVITY, g, USE_IMPLICIT_EULER);
        setup_mesh(args[0], q, qdot, V, F, alpha0, E, edge_target_angle, l0, edge_adjacent_vertices, k_axial, k_crease, EA, k_fold, k_facet, k_face, face_adjacent_edges);

        // Get path to activation profiles
        args[0].replace(args[0].find("crease_patterns"), std::string("crease_patterns").length(), "activation_profiles");
        args[0].replace(args[0].find(".fold"), 5, ".json");
        setup_dynamic_target_angles(args[0], edge_target_angle);
    } else {
        // No arguments provided, using default 
        std::cout << "No arguments provided, using default" << std::endl;
        setup_simulation_params("../data/simulation_params/default-params.json", dt, vertexMass, EA, k_fold, k_facet, k_face, zeta, ENABLE_STRAIN_VISUALIZATION, STRAIN_TYPE, ENABLE_DYNAMIC_SIMULATION, ENABLE_GRAVITY, g, USE_IMPLICIT_EULER);
        setup_mesh("../data/crease_patterns/defaultsquare.fold", q, qdot, V, F, alpha0, E, edge_target_angle, l0, edge_adjacent_vertices, k_axial, k_crease, EA, k_fold, k_facet, k_face, face_adjacent_edges);
        setup_dynamic_target_angles("../data/activation_profiles/defaultsquare.json", edge_target_angle);
    }

    // Set up mass matrix
    make_mass_matrix(M, q, vertexMass);

    // Set up previous angle tracking in crease constraints
    setup_prev_angle(edge_target_angle.size());

    // Create viewer
    igl::opengl::glfw::Viewer viewer;
    viewer_ptr = &viewer;

    // Set mesh and options
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.core().is_animating = true; 
    viewer.data().double_sided = true;
    viewer.core().lighting_factor = 0.0f;

    // Start the simulation in a seperate thread
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    if (ENABLE_STRAIN_VISUALIZATION){
        // Start the energy status in another thread if needed
        std::thread strain_visualization(visualizeStrain);
        strain_visualization.detach();
    }

    bool first_frame = true;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& v) -> bool {

        // the first frame needs to maximise the window. comment this out to disable automatic fullscreen
        if (first_frame){
            glfwMaximizeWindow(v.window);
            first_frame = false;
        }


        // This forces a redraw every frame
        return false;
    };

    // Launch viewer
    viewer.launch();
    
    // Clean up
    viewer_ptr = nullptr;
    simulating = false;

    if (args.size() >= 1){
        writeAverageStrainDuringSimulation(args[0]);
    } else {
        writeAverageStrainDuringSimulation("defaultsquare");
    }
    
    
    return 0;
}