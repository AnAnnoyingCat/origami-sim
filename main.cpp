// C++ and libraries
#include <igl/opengl/glfw/Viewer.h>
#include <chrono>
#include <thread>
#include <Eigen/Dense>

// Project Functions
#include <setup.h>
#include <drawHelperFunctions.h>
#include <assemble_edge_forces.h>
#include <d2V_axial_dq2.h>
#include <dV_axial_dq.h>
#include <forward_Euler.h>
#include <V_axial.h>
#include <T_vert.h>
#include <linearly_implicit_euler.h>
#include <assemble_stiffness.h>
#include <make_mass_matrix.h>

//Simulation state
bool simulating = true;
Eigen::VectorXd q;          // Generalized Vertex Coordinates, size 3*x vector
Eigen::VectorXd qdot;       // Generalized Vertex Velocities
Eigen::VectorXd x0;         // Fixed Point Constraints

Eigen::MatrixXd V;          // Vertices of the CP, nx3 matrix
Eigen::MatrixXi F;          // Faces of the CP, mx3 matrix
Eigen::MatrixXi E;          // Edges of the CP (Springs), ex2 matrix
Eigen::VectorXd l0;         // Original length of the Edges, size e vector
Eigen::VectorXd k;          // Per edge stiffness constant
Eigen::SparseMatrix<double> M; // Mass matrix

double yM = 1.0;            // Normalized Young's modulus
double csa = 100.0;         // Scaled cross-sectional area (Both used to calculate per axis stiffness)
double EA = yM * csa;       // Axial stiffness parameter
double t = 0;               // Simulation Time
double dt = 0.001;          // Time Step
double baseline_k = 1e5;    // Stiffness
double vertexMass = 1;      // Per vertex mass. Currently constant

// Working memory for integrator
Eigen::SparseMatrix<double> tmp_stiffness;
Eigen::VectorXd tmp_force;

//Pointer to the viewer
igl::opengl::glfw::Viewer* viewer_ptr = nullptr;


/// @brief Prints to console the total energy of the system - should be constant with no new energy introduced
void print_energy_status(){
    while(simulating){
        double V_edge, T_vertex, KE, PE;
        KE = 0;
        PE = 0;
        
        // Calculate per vertex kinetic energy
        for (int p = 0; p < V.rows(); p++) {
            T_vert(T_vertex, qdot.segment<3>(3 * p), vertexMass);

            KE += T_vertex;
        }

        // Calculate per axis potential energy
        for (int p = 0; p < E.rows(); p++){
            int v0 = E(p, 0);
            int v1 = E(p, 1);

            V_axial(V_edge, q.segment<3>(3 * v0), q.segment<3>(3 * v1), l0(p), k(p));

            PE += V_edge;
        }

        std::cout << "=====================================================" << std::endl;
        std::cout << "Kinetic Energy of the system: " << KE << std::endl;
        std::cout << "Potential Energy of the system: " << PE << std::endl;
        std::cout << "Total Energy of the system: " << KE + PE << std::endl;

        // Don't spam the console too much
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

/// @brief Main simulation loop, running in a seperate thread
void simulate(){
    while (simulating){

        auto forces = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            // Set f to zero and then add all the forces to it
            f.resize(q.size());
            f.setZero();

            assemble_edge_forces(f, q, E, l0, k);
            
            // TODO: ADD MORE FORCES -> Creases and Faces
        };

        auto stiffness = [&](Eigen::SparseMatrix<double> &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot) { 
            assemble_stiffness(K, q, qdot, V, E, l0, k);
        };

        
        //forward_euler(q, qdot, dt, forces, tmp_force);

        linearly_implicit_euler(q, qdot, dt, M, forces, stiffness, tmp_force, tmp_stiffness);

        // update vertex positions from q and increment time
        updateV(V, q);

        // Update the viewer's mesh
        if (viewer_ptr){
            viewer_ptr->data().set_vertices(V);
            viewer_ptr->data().compute_normals();
        }

        // Next time step
        t += dt;

        //std::cout << "qdot: " << qdot << std::endl;
        //simulating = false;

        // Small delay to make the animation visible
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

int main(int argc, char *argv[])
{   
    // Call setup to set up all the meshes and variables (may be changed later)
    setup(q, qdot, x0, V, F, E, l0, k, EA);
    make_mass_matrix(M, q, vertexMass);
    q(0) = 0.5;
    q(1) = 0.5;
    updateV(V, q);

    // Create viewer
    igl::opengl::glfw::Viewer viewer;
    viewer_ptr = &viewer;

    // Set mesh and options
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    viewer.core().is_animating = true; 

    // Start the simulation in a seperate thread
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    // Start the energy status in another thread
    std::thread total_energy_thread(print_energy_status);
    total_energy_thread.detach();
    
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) -> bool {
        // This forces a redraw every frame
        return false;
    };

    // Launch viewer
    viewer.launch();
    
    // Clean up
    viewer_ptr = nullptr;
    simulating = false;
    
    return 0;
}