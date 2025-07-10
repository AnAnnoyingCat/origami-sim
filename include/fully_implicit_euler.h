#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

/**
 * @brief Fully implicit Euler time integration scheme. 
 * 
 * @tparam FORCE 			
 * @tparam STIFFNESS 		
 * @param q 				Vertex coordinates
 * @param qdot 				Vertex velocities
 * @param dt 				Timestep size
 * @param mass 				Sparse mass matrix
 * @param force 			Function to calculate forces. Takes as input the vector to store the result, q and qdot.
 * @param stiffness 		Function to calculate stiffness. Takes as input the vector to store the result, q and qdot.
 * @param tmp_force 		Scratch space to calculate forces in
 * @param tmp_stiffness 	Scratch space to calculate stiffness in
 */
template<typename FORCE, typename STIFFNESS>
inline void implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, const Eigen::SparseMatrix<double> &mass, FORCE &force, STIFFNESS &stiffness, Eigen::VectorXd &tmp_force, Eigen::SparseMatrix<double> &tmp_stiffness, SimulationData& simulationData, SimulationParams& simulationParams) {
    
    const int MAX_ITERATIONS = 20;
    const double TOLERANCE = 1e-12;
    
    Eigen::VectorXd q_prev = q;
    Eigen::VectorXd qdot_prev = qdot;

    // Reuse the same vectors
    Eigen::VectorXd G(q.size());
    Eigen::VectorXd dq(q.size());

    for (int i = 0; i < MAX_ITERATIONS; i++){
        qdot = (q - q_prev) / dt;

        // Get the currnent deformed positions 
        get_deformed_positions(simulationData, simulationParams);

        // Detect active collisions
        simulationData.collisions.build(simulationData.collision_mesh, simulationData.deformed_vertices, simulationParams.min_barrier_distance);
        // Set up friction specific collision
        simulationData.friction_collisions.build(simulationData.collision_mesh, simulationData.deformed_vertices, simulationData.collisions, simulationData.barrier_potential, simulationData.barrier_stiffness, simulationParams.mu);

        force(tmp_force, q, qdot, (i == 0));
        stiffness(tmp_stiffness, q, qdot);

        G = mass * (q - q_prev - qdot_prev*dt) - tmp_force * dt * dt;
        
        if (G.norm() < TOLERANCE){
            break;
        }
        
        Eigen::SparseMatrix<double> J = mass - dt * dt * tmp_stiffness;

        // solve J * dq = -G
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(J);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Matrix decomposition failed");
        }
        dq = solver.solve(-G);


        // Update q
        q += dq;

    }

    // Update velocity
    qdot = (q - q_prev) / dt;
}