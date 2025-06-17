#include <Eigen/Dense>
#include <strain_calculations.h>
#include <map>
#include <trig_helper_functions.h>
#include <iostream>
#include <optional>

struct FoldInstruction {
    Interval time;    
    int fold_number;      
    double start_angle; 
    double end_angle;
    std::optional<std::string> mode; 
};

// Given a time t, updates edge_target_angle to that t
void calculateDynamicTargetAngle(Eigen::VectorXd &edge_target_angle, double t, Eigen::VectorXd q, Eigen::MatrixXi edge_adjacent_vertices);

/// @brief Passes all the time dependant target angle information to the function
void setupFoldTimeline(std::vector<FoldInstruction> in, double numEdges);

