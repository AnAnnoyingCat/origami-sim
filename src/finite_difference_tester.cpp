#include <iostream>
#include <Eigen/Dense>
#include <axial_constraints.h>
#include <axial_hessian.h>
#include <crease_constraints.h>
#include <crease_hessian.h>
#include <finite_difference_tester.h>
#include <assemble_damping_forces.h>
#include <assemble_damping_stiffness.h>
#include <assemble_barrier_forces.h>
#include <assemble_barrier_stiffness.h>

void test_axial_hessian(){
    double epsilon = 1e-6;
    Eigen::Vector3d q0(0, 0, 0);
    Eigen::Vector3d q1(1.25, 0, 0); // Deformed from rest length
    double l0 = 1;
    double stiffness = 20 / l0;

    
    Eigen::Matrix<double, 6, 1> q;
    q << q0, q1;

    // Storage for analytical derivative
    Eigen::Matrix<double, 6, 6> df_analytical;
    dF_axial(df_analytical, q0, q1, l0, stiffness);

    // Storage for numerical derivative
    Eigen::Matrix<double, 6, 6> df_numerical;
    df_numerical.setZero();

    // Central difference approximation
    for (int i = 0; i < 6; ++i) {
        Eigen::Matrix<double, 6, 1> q_perturbed_p = q;
        Eigen::Matrix<double, 6, 1> q_perturbed_n = q;

        q_perturbed_p(i) += epsilon;
        q_perturbed_n(i) -= epsilon;

        Eigen::Vector3d q0_p = q_perturbed_p.head<3>();
        Eigen::Vector3d q1_p = q_perturbed_p.tail<3>();
        Eigen::Vector3d q0_n = q_perturbed_n.head<3>();
        Eigen::Vector3d q1_n = q_perturbed_n.tail<3>();

        Eigen::Matrix<double, 6, 1> f_p, f_n;
        F_axial(f_p, q0_p, q1_p, l0, stiffness);
        F_axial(f_n, q0_n, q1_n, l0, stiffness);

        // Central difference: (f(x+e) - f(x-e)) / (2e)
        df_numerical.col(i) = (f_p - f_n) / (2 * epsilon);
    }

    // Compute and print the difference
    Eigen::Matrix<double, 6, 6> diff = df_numerical - df_analytical;
    double error = diff.norm();
	std::cout << "============== Running axial hessian finite difference test ============\n" << std::endl;
    std::cout << "Analytical derivative:\n" << df_analytical << std::endl;
    std::cout << "Numerical derivative:\n" << df_numerical << std::endl;
    std::cout << "Difference (numerical - analytical):\n" << diff << std::endl;
    std::cout << "Frobenius norm of difference: " << error << std::endl;
	std::cout << "========================================================================\n" << std::endl;
}

void test_crease_hessian(){
    double epsilon = 1e-6;
    Eigen::Vector3d q1(0, 0, 0);
    Eigen::Vector3d q2(1, 1, 0);
	Eigen::Vector3d q3(0, 1, 0);
	Eigen::Vector3d q4(1, 0, 0);
    double k_crease = 1.0;
    double theta_target = 180.0 * M_PI / 180; // 90°

    
    Eigen::Matrix<double, 12, 1> q;
    q << q1, q2, q3, q4;
	
    // Storage for analytical derivative
    Eigen::Matrix<double, 12, 12> df_analytical;
	df_analytical.setZero();
    dF_crease(df_analytical, q1, q2, q3, q4, k_crease, theta_target);

    // Storage for numerical derivative
    Eigen::Matrix<double, 12, 12> df_numerical;
    df_numerical.setZero();

    // Central difference approximation
    for (int i = 0; i < 12; i++) {
        Eigen::Matrix<double, 12, 1> q_perturbed_p = q;
        Eigen::Matrix<double, 12, 1> q_perturbed_n = q;

        q_perturbed_p(i) += epsilon;
        q_perturbed_n(i) -= epsilon;

        Eigen::Vector3d q1_p = q_perturbed_p.segment<3>(0);
        Eigen::Vector3d q2_p = q_perturbed_p.segment<3>(3);
		Eigen::Vector3d q3_p = q_perturbed_p.segment<3>(6);
		Eigen::Vector3d q4_p = q_perturbed_p.segment<3>(9);

        Eigen::Vector3d q1_n = q_perturbed_n.segment<3>(0);
        Eigen::Vector3d q2_n = q_perturbed_n.segment<3>(3);
		Eigen::Vector3d q3_n = q_perturbed_n.segment<3>(6);
		Eigen::Vector3d q4_n = q_perturbed_n.segment<3>(9);

        Eigen::Matrix<double, 12, 1> f_p, f_n;
		setup_prev_angle(2);
		
		F_crease(f_p, q1_p, q2_p, q3_p, q4_p, k_crease, theta_target, 0);
		F_crease(f_n, q1_n, q2_n, q3_n, q4_n, k_crease, theta_target, 1);

        // Central difference: (f(x+e) - f(x-e)) / (2e)
        df_numerical.col(i) = (f_p - f_n) / (2 * epsilon);
		
    }

    // Compute and print the difference
    Eigen::Matrix<double, 12, 12> diff = df_numerical - df_analytical;
    double error = diff.norm();

	std::cout << "============== Running crease hessian finite difference test ============\n" << std::endl;
    std::cout << "Analytical derivative:\n" << df_analytical << std::endl;
    std::cout << "Numerical derivative:\n" << df_numerical << std::endl;
    std::cout << "Difference (numerical - analytical):\n" << diff << std::endl;
    std::cout << "Frobenius norm of difference: " << error << std::endl;
	std::cout << "========================================================================\n" << std::endl;
}

void test_crease_force(){
    double epsilon = 1e-6;
    Eigen::Vector3d q1(0, 0, 0);
    Eigen::Vector3d q2(1, 1, 0);
    Eigen::Vector3d q3(0, 1, 0);
    Eigen::Vector3d q4(1, 0, 0);
    double k_crease = 1.0;
    double theta_target = 180.0 * M_PI / 180; // 90°

    Eigen::Matrix<double, 12, 1> q;
    q << q1, q2, q3, q4;

    // Storage for analytical derivative
    Eigen::Matrix<double, 12, 1> F_crease_by_me;
    F_crease_by_me.setZero();
    setup_prev_angle(1);
    F_crease(F_crease_by_me, q1, q2, q3, q4, k_crease, theta_target, 0);

    // Storage for numerical derivative
    Eigen::Matrix<double, 12, 1> F_crease_maple_res;
    F_crease_maple_res.setZero();
    F_crease_maple(F_crease_maple_res, q1, q2, q3, q4, k_crease, theta_target);

    // Compute and print the difference
    Eigen::Matrix<double, 12, 1> diff = F_crease_by_me - F_crease_maple_res;
    double error = diff.norm();

    std::cout << "============== Testing crease force and comparing maple ============\n" << std::endl;
    std::cout << "My force:\n" << F_crease_by_me << std::endl;
    std::cout << "Maple force:\n" << F_crease_maple_res << std::endl;
    std::cout << "Difference (mine - maple):\n" << diff << std::endl;
    std::cout << "Frobenius norm of difference: " << error << std::endl;
    std::cout << "========================================================================\n" << std::endl;
}

void test_damping_stiffness(){
    double epsilon = 1e-6;
    Eigen::Vector3d qdot1(0, 0, 0);
    Eigen::Vector3d qdot2(1, 1, 0);

    Eigen::Matrix<double, 6,1> qdot;
    qdot << qdot1, qdot2;

    Eigen::Matrix<int,1,2> E;
    E << 0,1;

    Eigen::VectorXd k_axial;
    k_axial.resize(1);
    k_axial << 1.0;
    double zeta = 4.5;

    

    // Storage for analytical derivative
    Eigen::SparseMatrix<double> stiffness_by_me;
    stiffness_by_me.setZero();
    stiffness_by_me.resize(6,6);
    assemble_damping_stiffness(stiffness_by_me, qdot, E, k_axial, zeta);

    // Storage for numerical derivative
    Eigen::Matrix<double, 6, 6> stiffness_numerical;
    stiffness_numerical.setZero();

    // Central difference approximation
    for (int i = 0; i < 6; i++) {
        Eigen::Matrix<double, 6, 1> qdot_perturbed_p = qdot;
        Eigen::Matrix<double, 6, 1> qdot_perturbed_n = qdot;

        qdot_perturbed_p(i) += epsilon;
        qdot_perturbed_n(i) -= epsilon;

        Eigen::VectorXd f_p, f_n;
        f_p.resize(6);
        f_n.resize(6);
        f_p.setZero();
        f_n.setZero();
		
		assemble_damping_forces(f_p, qdot_perturbed_p, E, k_axial, zeta);
        assemble_damping_forces(f_n, qdot_perturbed_n, E, k_axial, zeta);

        // Central difference: (f(x+e) - f(x-e)) / (2e)
        stiffness_numerical.col(i) = (f_p - f_n) / (2 * epsilon);

        // Debug prints
        std::cout << "Perturbation index: " << i << std::endl;
        std::cout << "qdot_perturbed_p:\n" << qdot_perturbed_p.transpose() << std::endl;
        std::cout << "qdot_perturbed_n:\n" << qdot_perturbed_n.transpose() << std::endl;
        std::cout << "f(q+ε): " << f_p.transpose() << std::endl;
        std::cout << "f(q-ε): " << f_n.transpose() << std::endl;
        std::cout << "Column " << i << " of numerical stiffness:\n" << stiffness_numerical.col(i).transpose() << "\n\n";
		
    }


    // Compute and print the difference
    Eigen::Matrix<double, 6, 6> diff = stiffness_by_me - stiffness_numerical;
    double error = diff.norm();

    std::cout << "============== Testing crease force and comparing maple ============\n" << std::endl;
    std::cout << "My force:\n" << stiffness_by_me << std::endl;
    std::cout << "Numerical force:\n" << stiffness_numerical << std::endl;
    std::cout << "Difference (mine - Numerical):\n" << diff << std::endl;
    std::cout << "Frobenius norm of difference: " << error << std::endl;
    std::cout << "========================================================================\n" << std::endl;
}

// void test_barrier_stiffness(){
//     double epsilon = 1e-5;
//     double mindistance = 1.0;
//     double k_barrier = 1.0;
//     Eigen::Vector3d q(0, 0, 0.2);

//     Eigen::VectorXd k_axial;
//     SimulationParams params;
//     params.k_barrier = 1.0;

//     // Storage for analytical derivative
//     Eigen::SparseMatrix<double> stiffness_by_me;
//     stiffness_by_me.setZero();
//     stiffness_by_me.resize(3,3);
//     assemble_ground_barrier_stiffness(stiffness_by_me, q, mindistance, k_barrier);

//     // Storage for numerical derivative
//     Eigen::Matrix<double, 3, 3> stiffness_numerical;
//     stiffness_numerical.setZero();

//     // Central difference approximation
//     for (int i = 0; i < 3; i++) {
//         Eigen::Vector3d q_perturbed_p = q;
//         Eigen::Vector3d q_perturbed_n = q;

//         q_perturbed_p(i) += epsilon;
//         q_perturbed_n(i) -= epsilon;

//         Eigen::VectorXd f_p, f_n;
//         f_p.resize(3);
//         f_n.resize(3);
//         f_p.setZero();
//         f_n.setZero();
		
// 		assemble_ground_barrier_forces(f_p, q_perturbed_p, mindistance, params);
//         assemble_ground_barrier_forces(f_n, q_perturbed_n, mindistance, params);

//         // Central difference: (f(x+e) - f(x-e)) / (2e)
//         stiffness_numerical.col(i) = (f_p - f_n) / (2 * epsilon);

//         // Debug prints
//         std::cout << "Perturbation index: " << i << std::endl;
//         std::cout << "qdot_perturbed_p:\n" << q_perturbed_p.transpose() << std::endl;
//         std::cout << "qdot_perturbed_n:\n" << q_perturbed_n.transpose() << std::endl;
//         std::cout << "f(q+ε): " << f_p.transpose() << std::endl;
//         std::cout << "f(q-ε): " << f_n.transpose() << std::endl;
//         std::cout << "Column " << i << " of numerical stiffness:\n" << stiffness_numerical.col(i).transpose() << "\n\n";
		
//     }


//     // Compute and print the difference
//     Eigen::Matrix<double, 3, 3> diff = stiffness_by_me - stiffness_numerical;
//     double error = diff.norm();

//     std::cout << "============== Testing barrier force and stiffness =====================\n" << std::endl;
//     std::cout << "My force:\n" << stiffness_by_me << std::endl;
//     std::cout << "Numerical force:\n" << stiffness_numerical << std::endl;
//     std::cout << "Difference (mine - Numerical):\n" << diff << std::endl;
//     std::cout << "Frobenius norm of difference: " << error << std::endl;
//     std::cout << "========================================================================\n" << std::endl;
// }