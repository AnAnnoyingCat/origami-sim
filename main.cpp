// C++ and libraries
#include <igl/opengl/glfw/Viewer.h>
#include <chrono>
#include <thread>
#include <Eigen/Dense>

// Project Functions
#include <setup.h>
#include <drawHelperFunctions.h>

//Simulation state
bool simulating = true;
Eigen::VectorXd q;          // Generalized Vertex Coordinates, size 3*x vector
Eigen::VectorXd qdot;       // Generalized Vertex Velocities
Eigen::VectorXd x0;         // Fixed Point Constraints

Eigen::MatrixXd V;          // Vertices of the CP, nx3 matrix
Eigen::MatrixXi F;          // Faces of the CP, mx3 matrix
Eigen::MatrixXi E;          // Edges of the CP (Springs), ex2 matrix
Eigen::VectorXd l0;         // Original length of the Edges, size e vector

double t = 0;               // Simulation Time
double dt = 0.005;          // Time Step
double k = 1e5;             // Stiffness

//Pointer to the viewer
igl::opengl::glfw::Viewer* viewer_ptr = nullptr;

void simulate(){
    while (simulating){
        // do the simulation 

        // Calculate Forces

        // Calculate Stiffness

        // Do time integration

        // update vertex positions from q and increment time
        updateV(V, q);

        // Update the viewer's mesh
        if (viewer_ptr){
            viewer_ptr->data().set_vertices(V);
            viewer_ptr->data().compute_normals();
        }
        t += dt;

        // Small delay to make the animation visible
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

int main(int argc, char *argv[])
{   
    // Call setup to set up all the meshes and variables (may be changed later)
    setup(q, qdot, x0, V, F, E, l0);

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