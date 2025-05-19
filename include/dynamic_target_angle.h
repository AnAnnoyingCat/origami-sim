#include <Eigen/Dense>
#include <strain_calculations.h>
#include <map>

struct FoldInstruction {
    Interval time;    
    int fold_number;      
    double start_angle; 
    double end_angle;
};

// Given a time t, updates edge_target_angle to that t
void calculateDynamicTargetAngle(Eigen::VectorXd& edge_target_angle, double t);

/// @brief Passes all the time dependant target angle information to the function
void setupTimeDependantTargetAngle(std::vector<FoldInstruction> in, double numEdges);

