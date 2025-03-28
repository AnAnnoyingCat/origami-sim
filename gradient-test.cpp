// C++ includes
#include <iostream>

// autodiff include
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include <Eigen/Dense>
using namespace autodiff;

void psi_neo_hookean(double &psi, Eigen::Ref<const Eigen::Matrix3d> F, double C, double D) {
    double Fdet = F.determinant();
    psi = C * (std::pow(Fdet, -2.0/3.0) * ((F.transpose() * F).trace()) - 3) + D * std::pow((Fdet - 1), 2);
}

// Custom determinant for 3x3 matrix with AutoDiff variables
var det3x3(const MatrixXvar& F) {
    return F(0,0)*(F(1,1)*F(2,2) - F(1,2)*F(2,1)) 
         - F(0,1)*(F(1,0)*F(2,2) - F(1,2)*F(2,0)) 
         + F(0,2)*(F(1,0)*F(2,1) - F(1,1)*F(2,0));
}

// real psi_neo_hookean_autodiff(const ArrayXreal& F, const real& C, const real& D){
	// // Reshape the 1D array F into a 3x3 matrix
	// Eigen::Map<const MatrixXreal> F_matrix(F.data(), 3, 3);

	// // Compute the determinant of F
	// real Fdet = F_matrix.determinant();

	// // Compute the Neo-Hookean energy
	// real psi = C * (pow(Fdet, -2.0/3.0) * ((F_matrix.transpose() * F_matrix).trace()) - 3) + D * pow(Fdet - 1, 2);

	// return psi;
// }

struct Params
{
    var C;
    var D;
};

var psi_neo_hookean_autodiff(const ArrayXvar& F, const Params& params){
	// Reshape the 1D array F into a 3x3 matrix
	Eigen::Map<const MatrixXvar> F_matrix(F.data(), 3, 3);
	// Compute the determinant of F
	var Fdet = det3x3(F_matrix);
	
	// Compute the Neo-Hookean energy
	var psi = params.C * (pow(Fdet, -2.0/3.0) * ((F_matrix.transpose() * F_matrix).trace()) - 3) + params.D * pow(Fdet - 1, 2);

	return psi;
}


int main()
{
	using Eigen::VectorXd;
	using Eigen::MatrixXd;

	// Baseline no autodiff
	Eigen::Matrix3d F;
	F << 1, 0, 0,
		0, 2, 0,
		0, 0, 1;
	double C = 1.2;
	double D = 1.5;
	double psi;
	psi_neo_hookean(psi, F, C, D);
	std::cout << "psi_neo_hookean without autodiff = " << psi << std::endl;
    
	// Base function autodiff
	ArrayXvar F_autodiff(9);
	F_autodiff << 1, 0, 0, 0, 2, 0, 0, 0, 1;
	Params params;
	params.C = 1.2;
	params.D = 1.5;

	std::cout << "psi_neo_hookean with autodiff = " << psi_neo_hookean_autodiff(F_autodiff, params) << std::endl;

	// gradient w/ autodiff
	var y = psi_neo_hookean_autodiff(F_autodiff, params);
	VectorXd gf = gradient(y, F_autodiff);

    std::cout << "gradient is = \n" << gf << std::endl;   

		// gradient with finite differences (1e-10)
		double ε = 5e-7;
		// F00 + ε
		Eigen::Matrix3d F00;
		F00 << 1 + ε, 0, 0,
			0, 2, 0,
			0, 0, 1;
		double dpsidf00;
		psi_neo_hookean(dpsidf00, F00, C, D);
		dpsidf00 = (dpsidf00 - psi) / ε;
	
		std::cout << "Derivative wrt f00 = " << dpsidf00 << std::endl;
		
		// F11 + ε
		Eigen::Matrix3d F11;
		F11 << 1, 0, 0,
			0, 2 + ε, 0,
			0, 0, 1;
		double dpsidf11;
		psi_neo_hookean(dpsidf11, F11, C, D);
		dpsidf11 = (dpsidf11 - psi) / ε;
		std::cout << "Derivative wrt f11 = " << dpsidf11 << std::endl;
		// F22 + ε
		Eigen::Matrix3d F22;
		F11 << 1, 0, 0,
			0, 2, 0,
			0, 0, 1 + ε;	

	// hessian w/ autodiff
	var u = psi_neo_hookean_autodiff(F_autodiff, params);
	//gf exists already
	Eigen::MatrixXd H = hessian(u, F_autodiff, gf);

	std::cout << "gradient (from hessian) is = \n" << gf << std::endl;   
	std::cout << "hessian is = \n" << H << std::endl;  
	
	Eigen::MatrixXd H_fd(9, 9);

	for (int i = 0; i < 9; i++) {
		for (int j = 0; j < 9; j++) {
			// Create perturbations in the i-th and j-th components
			ArrayXvar F_pp = F_autodiff;
			ArrayXvar F_pm = F_autodiff;
			ArrayXvar F_mp = F_autodiff;
			ArrayXvar F_mm = F_autodiff;
	
			F_pp(i) += ε;
			F_pp(j) += ε;
	
			F_pm(i) += ε;
			F_pm(j) -= ε;
	
			F_mp(i) -= ε;
			F_mp(j) += ε;
	
			F_mm(i) -= ε;
			F_mm(j) -= ε;
	
			// Evaluate psi at the perturbed points
			var psi_pp = psi_neo_hookean_autodiff(F_pp, params);
			var psi_pm = psi_neo_hookean_autodiff(F_pm, params);
			var psi_mp = psi_neo_hookean_autodiff(F_mp, params);
			var psi_mm = psi_neo_hookean_autodiff(F_mm, params);
	
			// Central difference formula
			H_fd(i, j) = (val(psi_pp) - val(psi_pm) - val(psi_mp) + val(psi_mm)) / (4.0 * ε * ε);
		}
	}
	
	std::cout << "Finite difference Hessian is = \n" << H_fd << std::endl;

}
