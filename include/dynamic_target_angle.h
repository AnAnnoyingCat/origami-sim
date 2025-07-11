#pragma once
#include <Eigen/Dense>
#include <strain_calculations.h>
#include <map>
#include <trig_helper_functions.h>
#include <iostream>
#include <optional>
#include <parameters.h>

struct FoldInstruction {
    Interval time;    
    int fold_number;      
    double start_angle; 
    double end_angle;
    std::optional<std::string> mode; 
};

// Given a time t, updates edge_target_angle to that t
void calculateDynamicTargetAngle(SimulationData& simulationData, SimulationParams& simulationParams);

/// @brief Passes all the time dependant target angle information to the function
void setupFoldTimeline(std::vector<FoldInstruction> in, double numEdges);

