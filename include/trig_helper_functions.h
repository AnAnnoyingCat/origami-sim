#include <Eigen/Dense>

void getAngle(double& alpha, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2);
void geth(double& h, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2);
void getNormal(Eigen::Vector3d& n, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, Eigen::Ref<const Eigen::Vector3d> q2);
void cot(double& alpha);