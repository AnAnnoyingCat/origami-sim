#include <Eigen/Dense>
#pragma once

/// @brief      Helper function to calculate per face strain based on internal angles. Strain scales linearly from 0% to 100% angle deformation (e.g. angle of 20 deg deformed to angle of 40 deg)
/// @param C    Color matrix for each face
void calculateFaceAngleStrain(Eigen::MatrixXd& C, Eigen::MatrixXi& F, Eigen::VectorXd& q, Eigen::MatrixXd& alpha0);

/// @brief      Helper function to calculate per axis strain, visualized by recoloring the faces of the mesh. TODO this isn't implemented yet
/// @param C    Color matrix for each face
void calculateAxialDeformationStrain(Eigen::MatrixXd& C, Eigen::MatrixXi& F, Eigen::MatrixXi E, Eigen::VectorXd& q, Eigen::VectorXd& l0, Eigen::MatrixXi& face_adjacent_edges);

/// @brief      Helper function to save the average strain of the entire execution to a file
void writeAverageStrainDuringSimulation();

/// @brief Helper class to describe an interval
class Interval {
    public:
    double begin;
    double end;

    Interval(double begin, double end) : begin(begin), end(end) {}

    // Define operator< for map key comparison
    bool operator<(const Interval& other) const {
        if (begin != other.begin)
            return begin < other.begin;
        return end < other.end;
    }

    bool contains(double query) const{
        return begin <= query && query < end;
    }
};