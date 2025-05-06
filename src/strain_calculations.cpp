#include <Eigen/Dense>
#include <strain_calculations.h>
#include <trig_helper_functions.h>

// Eigen::RowVector3d mesh_color(0.18, 0.78, 0.87);
Eigen::RowVector3d mesh_color_0(1.0, 1.0, 1.0);		// White 	-> 0  %
Eigen::RowVector3d mesh_color_1(0.0, 0.0, 1.0);		// Blue		-> 20 %
Eigen::RowVector3d mesh_color_2(0.0, 1.0, 0.0);		// Green 	-> 40 %
Eigen::RowVector3d mesh_color_3(1.0, 1.0, 0.0);		// Yellow	-> 60 %
Eigen::RowVector3d mesh_color_4(1.0, 0.0, 0.0);		// Red		-> 80 %
Eigen::RowVector3d mesh_color_5(0.7, 0.0, 1.0);		// Purple	-> 100%

void linearlyInterpretColors(Eigen::RowVector3d &C, Eigen::RowVector3d c1, Eigen::RowVector3d c2, double fraction){
	C = (1 - fraction) * c1 + fraction * c2;
}

void calculateFaceAngleStrain(Eigen::MatrixXd& C, Eigen::MatrixXi& F, Eigen::VectorXd& q, Eigen::MatrixXd& alpha0){
	C.resize(F.rows(), 3);

    for (int currFace = 0; currFace < F.rows(); currFace++){
        double currAlpha0, currAlpha1, currAlpha2;
        Eigen::Vector3d q0, q1, q2;
        q0 = q.segment<3>(3 * F(currFace, 0));
        q1 = q.segment<3>(3 * F(currFace, 1));
        q2 = q.segment<3>(3 * F(currFace, 2));
        
        getAngle(currAlpha0, q0, q1, q2);
        getAngle(currAlpha1, q1, q2, q0);
        getAngle(currAlpha2, q2, q0, q1);

        // Average strain of the angles
        //double strain = (1.0/3.0) * (std::abs(currAlpha0 - alpha0(currFace, 0)) / alpha0(currFace, 0) + std::abs(currAlpha1 - alpha0(currFace, 1)) / alpha0(currFace, 1) + std::abs(currAlpha2 - alpha0(currFace, 2)) / alpha0(currFace, 2));
        
        // Largest strain of the angles
        double strain = std::max(std::abs(currAlpha0 - alpha0(currFace, 0)) / alpha0(currFace, 0), std::max(std::abs(currAlpha1 - alpha0(currFace, 1)) / alpha0(currFace, 1), std::abs(currAlpha2 - alpha0(currFace, 2)) / alpha0(currFace, 2)));

        strain = std::min(strain, 1.0); // Safety check, unnecessary i think
		
        Eigen::RowVector3d faceColor;
		
		if (strain < 0.05) {
			// *5 is to scale the strain from a 0-20% range to a 0-100% range
			linearlyInterpretColors(faceColor, mesh_color_0, mesh_color_1, strain * 20);
		} else if (strain < 0.10) {
			linearlyInterpretColors(faceColor, mesh_color_1, mesh_color_2, (strain - 0.05) * 20);
		} else if (strain < 0.15){
			linearlyInterpretColors(faceColor, mesh_color_2, mesh_color_3, (strain - 0.10) * 20);
		} else if (strain < 0.20){
			linearlyInterpretColors(faceColor, mesh_color_3, mesh_color_4, (strain - 0.15) * 20);
		} else {
			linearlyInterpretColors(faceColor, mesh_color_4, mesh_color_5, (strain - 0.20) * 20);
		}

        C.row(currFace) = faceColor;
    }
}

/// @brief      Helper function to calculate per axis strain, visualized by recoloring the faces of the mesh. TODO this isn't implemented yet
/// @param C    Color matrix for each face
void calculateAxialDeformationStrain(Eigen::MatrixXd& C){

}