// C++ includes
#include <iostream>

// autodiff include
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include <Eigen/Dense>
using namespace autodiff;

void psi_neo_hookean(double &psi, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    double Fdet = F.determinant();
    psi = C * (std::pow(Fdet, -2.0/3.0) * ((F.transpose() * F).trace()) - 3) + D * std::pow((Fdet - 1), 2);
}

real psi_neo_hookean_autodiff(const ArrayXreal& F, const real& C, const real& D){
	// Reshape the 1D array F into a 3x3 matrix
	Eigen::Map<const MatrixXreal> F_matrix(F.data(), 3, 3);

	// Compute the determinant of F
	real Fdet = F_matrix.determinant();

	// Compute the Neo-Hookean energy
	real psi = C * (pow(Fdet, -2.0/3.0) * ((F_matrix.transpose() * F_matrix).trace()) - 3) + D * pow(Fdet - 1, 2);

	return psi;
}

int main()
{
	using Eigen::VectorXd;
	using Eigen::MatrixXd;

	Eigen::Matrix3d F;
	F << 1, 0, 0,
		0, 2, 0,
		0, 0, 1;
	double C = 1.2;
	double D = 1.5;
	double psi;
	psi_neo_hookean(psi, F, C, D);
	std::cout << "psi_neo_hookean without autodiff = " << psi << std::endl;
    
	ArrayXreal F_autodiff(9);
	F_autodiff << 1, 0, 0, 0, 2, 0, 0, 0, 1;
	real C_autodiff = 1.2;
	real D_autodiff = 1.5;

	std::cout << "psi_neo_hookean with autodiff = " << psi_neo_hookean_autodiff(F_autodiff, C_autodiff	, D_autodiff) << std::endl;

	real u;
	VectorXd gf;
	gradient(psi_neo_hookean_autodiff, wrt(F_autodiff), at(F_autodiff, C_autodiff, D_autodiff), u, gf);

	std::cout << "u = " << u << std::endl;       // print the evaluated output u
    std::cout << "gx = \n" << gf << std::endl;   

	MatrixXd Jf;
	VectorXreal Q;
	jacobian(psi_neo_hookean_autodiff, wrt(F_autodiff), at(F_autodiff, C_autodiff, D_autodiff), Q, Jf);

	std::cout << "gf = \n" << gf << std::endl;     // print the evaluated output vector F
    std::cout << "Jf = \n" << Jf << std::endl;

}
