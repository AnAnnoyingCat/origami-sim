#include <Eigen/Dense>
#include <strain_calculations.h>
#include <trig_helper_functions.h>
#include <map>
#include <iostream>
#include <fstream>

std::vector<double> frameStrains;
std::map<Interval, Eigen::Matrix<double, 6, 1>> colorMap{
    {Interval(0.00, 0.17), (Eigen::Matrix<double, 6, 1>() << 14,14,120, 62,117,207).finished()}, 
    {Interval(0.17, 0.30), (Eigen::Matrix<double, 6, 1>() << 62,117,207, 91,190,243).finished()},
    {Interval(0.30, 0.43), (Eigen::Matrix<double, 6, 1>() << 91,190,243, 175,237,234).finished()},
    {Interval(0.43, 0.50), (Eigen::Matrix<double, 6, 1>() << 175,237,234, 229,241,196).finished()},
    {Interval(0.50, 0.59), (Eigen::Matrix<double, 6, 1>() << 229,241,196, 244,213,130).finished()},
    {Interval(0.59, 0.71), (Eigen::Matrix<double, 6, 1>() << 244,213,130, 237,158,80).finished()},
    {Interval(0.71, 0.85), (Eigen::Matrix<double, 6, 1>() << 237,158,80, 204,90,41).finished()},
    {Interval(0.85, 1.01), (Eigen::Matrix<double, 6, 1>() << 204,90,41, 150,20,30).finished()}
};


/// @brief      Helper function to quickly lerp the color based on a provided colorMap
/// @param C    Color gets returned here
/// @param x    The provided x position
void getColor(Eigen::RowVector3d &C, double x){
    // Clamp x to the interval 0, 1
    x = std::min(std::max(0.0, x), 1.0);
    //std::cout << "x is: " << x << std::endl;

    // Find which interval x is in and get the respective color map
    Interval i(0,0);
    Eigen::Matrix<double, 6, 1> lerpColors;
    for (const auto& [interval, color] : colorMap) {
        //std::cout << "current interval is: " << interval.begin << ", " << interval.end << ", contains returns: " << interval.contains(x) << std::endl;
        if (interval.contains(x)) {
            i = interval;
            lerpColors = color;
            break;
        }
    }

    // Lineraly interpolate between the two colors
    double fraction = (x - i.begin) / (i.end - i.begin);
	C = (1 - fraction) * (lerpColors.segment<3>(0) / 255.0) + fraction * (lerpColors.segment<3>(3) / 255.0);
    //std::cout << "Fraction is: " << fraction << ", and final color is: " << C << std::endl;
}

void calculateFaceAngleStrain(Eigen::MatrixXd& C, Eigen::MatrixXi& F, Eigen::VectorXd& q, Eigen::MatrixXd& alpha0){
	C.resize(F.rows(), 3);

    double frameStrain = 0.0;
    for (int currFace = 0; currFace < F.rows(); currFace++){
        double currAlpha0, currAlpha1, currAlpha2;
        Eigen::Vector3d q0, q1, q2;
        q0 = q.segment<3>(3 * F(currFace, 0));
        q1 = q.segment<3>(3 * F(currFace, 1));
        q2 = q.segment<3>(3 * F(currFace, 2));
        
        getAngle(currAlpha0, q0, q1, q2);
        getAngle(currAlpha1, q1, q2, q0);
        getAngle(currAlpha2, q2, q0, q1);

        // Largest strain of the angles
        double strain = std::max(std::abs(currAlpha0 - alpha0(currFace, 0)) / alpha0(currFace, 0), std::max(std::abs(currAlpha1 - alpha0(currFace, 1)) / alpha0(currFace, 1), std::abs(currAlpha2 - alpha0(currFace, 2)) / alpha0(currFace, 2)));

        // By multiplying by 5 we effectively scale the strain from the interval 0 - 0.2 to the full interval 0 - 1
        strain = std::min(strain * 5, 1.0); // Safety check, unnecessary i think
		frameStrain += strain;
        Eigen::RowVector3d faceColor;

        getColor(faceColor, strain);
		
        C.row(currFace) = faceColor;
    }
    frameStrain /= F.rows();
    frameStrains.push_back(frameStrain);
}


void calculateAxialDeformationStrain(Eigen::MatrixXd& C, Eigen::MatrixXi& F, Eigen::MatrixXi E, Eigen::VectorXd& q, Eigen::VectorXd& l0, Eigen::MatrixXi& face_adjacent_edges){
    C.resize(F.rows(), 3);

    Eigen::VectorXd edgeStrains;
    edgeStrains.resize(E.rows());

    // Pre-Calc strain for each edge so we don't need to calculate it twice
    for (int currEdge = 0; currEdge < E.rows(); currEdge++){
        Eigen::Vector3d q0, q1;
        q0 = q.segment<3>(3 * E(currEdge, 0));
        q1 = q.segment<3>(3 * E(currEdge, 1));
        double currLen = (q1 - q0).norm();
        edgeStrains(currEdge) = (currLen - l0(currEdge)) / l0(currEdge);
    }

    // Then go over every face and assign it the largest edge deformation strain 
    for (int currFace = 0; currFace < F.rows(); currFace++){
        // Get the largest strain of the current face acjacent edges
        double strain = std::max(std::abs(edgeStrains(face_adjacent_edges(currFace, 0))), std::max(std::abs(edgeStrains(face_adjacent_edges(currFace, 1))), std::abs(edgeStrains(face_adjacent_edges(currFace, 2)))));

        // Multiply by 5 to effectively stretch the 0 - 0.2 interval to the interval 0 - 1. 
        strain = std::min(strain * 5, 1.0);

        Eigen::RowVector3d faceColor;
        
        getColor(faceColor, strain);
        C.row(currFace) = faceColor;
    }
    frameStrains.push_back(edgeStrains.cwiseAbs().mean());
}

void writeAverageStrainDuringSimulation(){
    double averageStrain = 0.0;
    for (double s : frameStrains) {
        averageStrain += s;
    }
    averageStrain /= frameStrains.size();

    std::ofstream outFile("../average_strain.txt");
    if (outFile.is_open()) {
        outFile << "Average strain over simulation: " << averageStrain << std::endl;
        outFile.close();
    } else {
        std::cerr << "Failed to open output file.\n";
    }
}