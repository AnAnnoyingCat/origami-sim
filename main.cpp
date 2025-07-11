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
#include <parameters.h>
#include <finite_difference_tester.h>
#include <assemble_barrier_forces.h>
#include <assemble_barrier_stiffness.h>
#include <assemble_friction_forces.h>
#include <assemble_friction_stiffness.h>
#include <IPC-helperfunctions.h>
#include <ipc/ipc.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/potentials/friction_potential.hpp>


// This keeps track of all data of the simulation
SimulationData simulationData;

// This keeps track of the parameters regarding the simulation
SimulationParams simulationParams;

// Working memory for the integrator
Eigen::SparseMatrix<double> tmp_stiffness;
Eigen::VectorXd tmp_force;

// Pointer to the viewer for update reasons
igl::opengl::glfw::Viewer* viewer_ptr = nullptr;
int simulation_mesh_id;

/// @brief If enabled in config, visualizes strain in the mesh using l0 or alpha0
void visualizeStrain(){
    while(simulationParams.simulating){
        if (viewer_ptr){
            Eigen::MatrixXd C;
            if (simulationParams.STRAIN_TYPE == "face"){
                calculateFaceAngleStrain(C, simulationData.F, simulationData.q, simulationData.alpha0);
            } else if (simulationParams.STRAIN_TYPE == "edge"){
                calculateAxialDeformationStrain(C, simulationData.F, simulationData.E, simulationData.q, simulationData.l0, simulationData.face_adjacent_edges);
            }
            viewer_ptr->data(simulation_mesh_id).set_colors(C);
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
    center(2) = simulationParams.spawn_height;
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
    while (simulationParams.simulating){
        
        // If simulation type is dynamic, recalculate the target angle for the curren frame
        if (simulationParams.ENABLE_DYNAMIC_SIMULATION){
            calculateDynamicTargetAngle(simulationData, simulationParams);
        }

        auto forces = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, bool first_time){
            // Set f to zero and then add all the forces to it
            f.resize(q.size());
            f.setZero();

            // Get the basic forces
            std::cout << std::fixed;
            std::cout << std::setprecision(3);
            assemble_edge_forces(f, q, simulationData.E, simulationData.l0, simulationData.k_axial);
            if (simulationParams.LOG_FORCES){
                std::cout << "Edge force: " << f(2);  
            }                                        
                            
            assemble_crease_forces(f, q, simulationData.edge_adjacent_vertices, simulationData.k_crease, simulationData.edge_target_angle); 
            if (simulationParams.LOG_FORCES){
                std::cout << " + Crease force: " << f(2);
            }     
            
            assemble_damping_forces(f, qdot, simulationData.E, simulationData.k_axial, simulationParams.zeta);
            if (simulationParams.LOG_FORCES){
                std::cout << " + damping force: " << f(2);
            }

            // If gravity is enabled, get that too
            if (simulationParams.ENABLE_GRAVITY){
                assemble_gravity_forces(f, simulationParams.g, simulationParams.vertexMass);
                if (simulationParams.LOG_FORCES){
                    std::cout << " + gravity force: " << f(2);
                }    
            }
            
            if (simulationParams.enable_barrier){
                assemble_barier_forces_IPC(f, simulationParams, simulationData, first_time);
                if (simulationParams.LOG_FORCES){
                    std::cout << " + Barrier force: " << f(2);
                }
            }

            if (simulationParams.enable_friction){
                assemble_friction_forces_IPC(f, simulationParams, simulationData);
                if (simulationParams.LOG_FORCES){
                    std::cout << " + friction force: " << f(2);
                }
            }

            if (simulationParams.LOG_FORCES){
                std::cout << std::endl;    
            }
        };

        auto stiffness = [&](Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            K.resize(q.size(), q.size());
            K.setZero();

            // Get all the basic stiffnessess
            assemble_edge_stiffness(K, q, simulationData.V, simulationData.E, simulationData.l0, simulationData.k_axial);
            assemble_crease_stiffness(K, q, simulationData.edge_adjacent_vertices, simulationData.k_crease, simulationData.edge_target_angle);
            assemble_damping_stiffness(K, qdot, simulationData.E, simulationData.k_axial, simulationParams.zeta);
            if (simulationParams.enable_barrier){
                assemble_barier_stiffness_IPC(K, simulationParams, simulationData);
            }
            if (simulationParams.enable_friction){
                assemble_friction_stiffness_IPC(K, simulationParams, simulationData);
            }
            
        };

        // Get the current collision mesh before caluculating all the forces
        make_collision_mesh(simulationData, simulationParams);

        // Time integration
        if (simulationParams.USE_IMPLICIT_EULER){
            implicit_euler(simulationData.q, simulationData.qdot, simulationParams.dt, simulationData.M, forces, stiffness, tmp_force, tmp_stiffness, simulationData, simulationParams);
        } else {
            forward_euler(simulationData.q, simulationData.qdot, simulationParams.dt, forces, tmp_force, simulationData, simulationParams);
        }
        

        // Make sure the floating mesh doesn't drift off with gravity disabled
        if (!simulationParams.ENABLE_GRAVITY){
            centerMesh();
        }

        // update vertex positions from q and increment time
        updateV(simulationData.V, simulationData.q);
        
        // Update the viewer's mesh
        if (viewer_ptr){
            viewer_ptr->data(simulation_mesh_id).set_vertices(simulationData.V);
            viewer_ptr->data(simulation_mesh_id).compute_normals();
        }

        // Next time step
        simulationData.t += simulationParams.dt;
        if (simulationParams.LOG_SIMULATION_TIME){
            std::cout << "t: " << simulationData.t << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{    

    // Read args into a vector
    std::vector<std::string> args;
    std::copy(argv + 1, argv + argc, std::back_inserter(args));

    // Set up parameters and read CP
    if (args.size() == 3){
        // Both CP and params provided
        setup_simulation_params(args[1], simulationParams);
        setup_mesh(args[0], simulationParams, simulationData);

        // Get path to activation profiles
        setup_dynamic_target_angles(args[2], simulationData.edge_target_angle);
    } else if (args.size() == 2){
        // Both CP and params provided, checking for actuation profiles at default path and name
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
    simulation_mesh_id = viewer.append_mesh();
    viewer.data(simulation_mesh_id).set_mesh(simulationData.V, simulationData.F);
    viewer.data(simulation_mesh_id).set_face_based(true);
    viewer.data(simulation_mesh_id).double_sided = true;
    viewer.core().is_animating = true; 
    viewer.core().lighting_factor = 0.0f;

    if (simulationParams.enable_floor){
        // Set up the floor mesh
        int floor_id = setup_floor(viewer, simulationData);

        // Adjust camera zoom
        viewer.core().camera_zoom = simulationParams.sim_zoom_level;
    }

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
        // This forces a redraw every frame
        return false;
    };

    // Launch viewer
    viewer.launch();
    
    // Clean up
    viewer_ptr = nullptr;
    simulationParams.simulating = false;

    if (args.size() >= 1){
        writeAverageStrainDuringSimulation(args[0]);
    } else {
        writeAverageStrainDuringSimulation("defaultsquare");
    }
    
    return 0;
}