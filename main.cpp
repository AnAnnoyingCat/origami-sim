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
#include <V_axial.h>
#include <T_vert.h>
#include <linearly_implicit_euler.h>
#include <assemble_stiffness.h>
#include <make_mass_matrix.h>
#include <assemble_damping_forces.h>

// Simulation state
bool simulating = true;
Eigen::VectorXd q;                  // Generalized Vertex Coordinates, size 3*x vector
Eigen::VectorXd qdot;               // Generalized Vertex Velocities
Eigen::SparseMatrix<double> M;      // Sparse Mass Matrix.
Eigen::SparseMatrix<double> P;      // Fixed Point Constraints
Eigen::VectorXd x0;                 // Fixed vertex indices

Eigen::MatrixXd V;          // Vertices of the CP, nx3 matrix
Eigen::MatrixXi F;          // Faces of the CP, mx3 matrix
Eigen::MatrixXi E;          // Edges of the CP (Springs), ex2 matrix
Eigen::VectorXd edge_theta; // Final Target fold angle for each edge
Eigen::VectorXd curr_theta; // Current Target fold angle for each edge (CURRENTLY UNUSED)
Eigen::VectorXd l0;         // Original length of the Edges, size e vector
Eigen::MatrixXd alpha0;     // Nominal angles in flat state of CP
Eigen::VectorXd k_axial;    // Per edge stiffness constant
Eigen::VectorXd k_crease;   // Per crease stiffness constant

std::vector<std::array<int, 4>> edge_adjacent_vertices; // For each edge, stores the four vertices making up the two triangles which meet at the edge. The order is: Right vertex, Left Vertex, Start Vertex, End Vertex. It is {-1, -1, -1, -1} for border edges

double t = 0;               // Simulation Time
double dt = 0.001;          // Time Step
double vertexMass = 1;      // Per vertex mass. Currently constant
double EA = 1.0 * 5e4;      // Axial stiffness parameter, used in calculating axial stiffness
double k_fold = 1e3;        // Stiffness for a mountain or valley crease (Should be much smaller than the axial stiffness)
double k_facet = 1e3;       // Stiffness for a facet crease
double k_face = 2e2;        // Stiffness for the face constraints
double zeta = 0.25;         // Parameter in the damping ratio from the paper

// Working memory for integrator
Eigen::SparseMatrix<double> tmp_stiffness;
Eigen::VectorXd tmp_force;

// Pointer to the viewer for update reasons
igl::opengl::glfw::Viewer* viewer_ptr = nullptr;

// Debug flags
const bool PRINT_FORCE_INFO = true;
const bool ANIMATE_ANGLE = true;

/// @brief Prints to console the total energy of the system - should be constant with no new energy introduced
void print_energy_status(){
    while(simulating){
        double V_edge, T_vertex, KE, PE;
        KE = 0;
        PE = 0;
        
        // Calculate per vertex kinetic energy
        for (int p = 0; p < V.rows(); p++) {
            T_vert(T_vertex, (P.transpose() * qdot).segment<3>(3 * p), vertexMass);

            KE += T_vertex;
        }

        // Calculate per axis potential energy
        for (int p = 0; p < E.rows(); p++){
            int v0 = E(p, 0);
            int v1 = E(p, 1);

            V_axial(V_edge, (P.transpose() * q + x0).segment<3>(3 * v0), (P.transpose() * q + x0).segment<3>(3 * v1), l0(p), k_axial(p));

            PE += V_edge;
        }

        std::cout << "=====================================================" << std::endl;
        std::cout << "Kinetic Energy of the system: " << KE << std::endl;
        std::cout << "Potential Energy of the system: " << PE << std::endl;
        std::cout << "Total Energy of the system: " << KE + PE << std::endl;

        // Don't spam the console too much
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }
}

/// @brief Helper function which updates V using q
void updateV(Eigen::MatrixXd& V, Eigen::VectorXd q){
	int n = V.rows();  
	for (int i = 0; i < n; i++){
		V.row(i) = q.segment<3>(3 * i);
	}
}

/// @brief Main simulation loop, running in a seperate thread
void simulate(){
    while (simulating){

        auto forces = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot){
            // Set f to zero and then add all the forces to it
            f.resize(q.size());
            f.setZero();
            assemble_edge_forces(f, P.transpose() * q + x0, E, l0, k_axial);
            assemble_crease_forces(f, P.transpose() * q + x0, edge_adjacent_vertices, k_crease, edge_theta);

            // TODO: assemble_face_forces()

            assemble_damping_forces(f, qdot, E, k_axial, zeta);
        };

        forward_euler(q, qdot, dt, forces, tmp_force);

        // update vertex positions from q and increment time
        updateV(V, P.transpose() * q + x0);
        

        // Update the viewer's mesh
        if (viewer_ptr){
            viewer_ptr->data().set_vertices(V);
            viewer_ptr->data().compute_normals();
        }

        // Next time step
        t += dt;

        // Small delay to make the animation visible
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

void animated_angle(){
    while (simulating){
        edge_theta(4) = -90.0 * M_PI / 180;
        std::this_thread::sleep_for(std::chrono::milliseconds(2500));
        edge_theta(4) = 90.0 * M_PI / 180;
        std::this_thread::sleep_for(std::chrono::milliseconds(2500));
    }
}

int main(int argc, char *argv[])
{   
    // Call setup to set up all the meshes and variables
    setup(q, qdot, x0, P, V, F, alpha0, E, edge_theta, l0, k_axial, k_crease, EA, k_fold, k_facet, edge_adjacent_vertices);

    // Set up mass matrix
    make_mass_matrix(M, q, vertexMass);
    M = P*M*P.transpose();

    // Introduce offset to check how forces react
    // q(0) = -1;
    // q(1) = -1;
    // updateV(V, P.transpose() * q + x0);
    

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

    if (PRINT_FORCE_INFO){
        // Start the energy status in another thread if needed
        std::thread total_energy_thread(print_energy_status);
        total_energy_thread.detach();
    }

    if (ANIMATE_ANGLE){
        // start the animation thread seperately. crude solution for now
        std::thread animation(animated_angle);
        animation.detach();
    }

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