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

// Set this to false to immediately stop all processes to do with simulation
bool simulating = true;              

// This keeps track of all data of the simulation
SimulationData simulationData;

// This keeps track of the parameters regarding the simulation
SimulationParams simulationParams;

// Working memory for the integrator
Eigen::SparseMatrix<double> tmp_stiffness;
Eigen::VectorXd tmp_force;

// Pointer to the viewer for update reasons
igl::opengl::glfw::Viewer* viewer_ptr = nullptr;

/// @brief If enabled in config, visualizes strain in the mesh using l0 or alpha0
void visualizeStrain(){
    while(simulating){
        if (viewer_ptr){
            Eigen::MatrixXd C;
            if (simulationParams.STRAIN_TYPE == "face"){
                calculateFaceAngleStrain(C, simulationData.F, simulationData.q, simulationData.alpha0);
            } else if (simulationParams.STRAIN_TYPE == "edge"){
                calculateAxialDeformationStrain(C, simulationData.F, simulationData.E, simulationData.q, simulationData.l0, simulationData.face_adjacent_edges);
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
    for (int currVertex = 0; currVertex < simulationData.V.rows(); currVertex++){
        center += simulationData.q.segment<3>(3 * currVertex);
    }
    center /= simulationData.V.rows();
    for (int currVertex = 0; currVertex < simulationData.V.rows(); currVertex++){
        simulationData.q.segment<3>(3 * currVertex) -= center;
    }
}

/// @brief Main simulation loop, running in a seperate thread
void simulate(){
    while (simulating){
        
        // If simulation type is dynamic, recalculate the target angle for the curren frame
        if (simulationParams.ENABLE_DYNAMIC_SIMULATION){
            calculateDynamicTargetAngle(simulationData.edge_target_angle, simulationData.t, simulationData.q, simulationData.edge_adjacent_vertices);
        }

        auto forces = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            // Set f to zero and then add all the forces to it
            f.resize(q.size());
            f.setZero();

            // Get the basic forces
            assemble_edge_forces(f, q, simulationData.E, simulationData.l0, simulationData.k_axial);                                         // Keep edge lengths constant
            assemble_damping_forces(f, qdot, simulationData.E, simulationData.k_axial, simulationParams.zeta);                                 // Apply viscous dampening
            assemble_crease_forces(f, q, simulationData.edge_adjacent_vertices, simulationData.k_crease, simulationData.edge_target_angle);  // Calculate crease driving force

            // If graivty is enabled, get that too
            if (simulationParams.ENABLE_GRAVITY){
                assemble_gravity_forces(f, simulationParams.g, simulationParams.vertexMass);
            }
        };

        auto stiffness = [&](Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            K.resize(q.size(), q.size());
            K.setZero();

            // Get all the basic stiffnessess
            assemble_edge_stiffness(K, q, simulationData.V, simulationData.E, simulationData.l0, simulationData.k_axial);
            assemble_crease_stiffness(K, q, simulationData.edge_adjacent_vertices, simulationData.k_crease, simulationData.edge_target_angle);
            assemble_damping_stiffness(K, qdot, simulationData.E, simulationData.k_axial, simulationParams.zeta);

            /* If gravity is enabled, and uses a more complicated than normal computation, get that too
            if (ENABLE_GRAVITY){
                assemble_gravity_stiffness(K);
            }
            */
            
        };

        // Time integration
        if (simulationParams.USE_IMPLICIT_EULER){
            implicit_euler(simulationData.q, simulationData.qdot, simulationParams.dt, simulationData.M, forces, stiffness, tmp_force, tmp_stiffness);
        } else {
            forward_euler(simulationData.q, simulationData.qdot, simulationParams.dt, forces, tmp_force);
        }
        

        // Make sure the floating mesh doesn't drift off with gravity disabled
        if (!simulationParams.ENABLE_GRAVITY){
            centerMesh();
        }

        // update vertex positions from q and increment time
        updateV(simulationData.V, simulationData.q);
        
        // Update the viewer's mesh
        if (viewer_ptr){
            viewer_ptr->data().set_vertices(simulationData.V);
            viewer_ptr->data().compute_normals();
        }

        // Next time step
        simulationData.t += simulationParams.dt;
        std::cout << "t: " << simulationData.t << std::endl;
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
        setup_simulation_params(args[1], simulationParams);
        setup_mesh(args[0], simulationParams, simulationData);

        // Get path to activation profiles
        args[0].replace(args[0].find("crease_patterns"), std::string("crease_patterns").length(), "activation_profiles");
        args[0].replace(args[0].find(".fold"), 5, ".json");
        setup_dynamic_target_angles(args[0], simulationData.edge_target_angle);

    } else if (args.size() == 1){
        // Only crease pattern provided
        std::cout << "Using default parameters" << std::endl;
        setup_simulation_params("../data/simulation_params/default-params.json", simulationParams);
        setup_mesh(args[0], simulationParams, simulationData);

        // Get path to activation profiles
        args[0].replace(args[0].find("crease_patterns"), std::string("crease_patterns").length(), "activation_profiles");
        args[0].replace(args[0].find(".fold"), 5, ".json");
        setup_dynamic_target_angles(args[0], simulationData.edge_target_angle);
    } else {
        // No arguments provided, using default 
        std::cout << "No arguments provided, using default" << std::endl;
        setup_simulation_params("../data/simulation_params/default-params.json", simulationParams);
        setup_mesh("../data/crease_patterns/defaultsquare.fold", simulationParams, simulationData);
        setup_dynamic_target_angles("../data/activation_profiles/defaultsquare.json", simulationData.edge_target_angle);
    }

    // Set up mass matrix
    make_mass_matrix(simulationData.M, simulationData.q, simulationParams.vertexMass);

    // Set up previous angle tracking in crease constraints
    setup_prev_angle(simulationData.edge_target_angle.size());

    // Create viewer
    igl::opengl::glfw::Viewer viewer;
    viewer_ptr = &viewer;

    // Set mesh and options
    viewer.data().set_mesh(simulationData.V, simulationData.F);
    viewer.data().set_face_based(true);
    viewer.core().is_animating = true; 
    viewer.data().double_sided = true;
    viewer.core().lighting_factor = 0.0f;

    // Start the simulation in a seperate thread
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    if (simulationParams.ENABLE_STRAIN_VISUALIZATION){
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