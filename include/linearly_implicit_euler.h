#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

// Input:
//  q - generalized coordinates for the mass-spring system
//  qdot - generalized velocity for the mass-spring system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(q, qdot) - a function that computes the force acting on the mass-spring system. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
// Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS>
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                                   const Eigen::SparseMatrix<double> &mass, FORCE &force, STIFFNESS &stiffness,
                                   Eigen::VectorXd &tmp_force, Eigen::SparseMatrix<double> &tmp_stiffness) {
    // Compute the stiffness Matrix
    stiffness(tmp_stiffness, q, qdot);

    // Compute the force vector
    force(tmp_force, q, qdot);

    // Form the linear system: (M - dt^2 * K) * qdot_{t+1} = M * qdot_t + dt * f(q_t)
    Eigen::SparseMatrix<double> A = mass - dt * dt * tmp_stiffness;
    Eigen::VectorXd b = mass * qdot + dt * tmp_force;

    // Solve the linear system for qdot_{t+1}
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(A);
    Eigen::VectorXd qdot_t1 = solver.solve(b);

    // Update q and qdot
    qdot = qdot_t1; // Update velocity
    q += dt * qdot; // Update position
}