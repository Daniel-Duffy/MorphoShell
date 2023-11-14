
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Node_Class.hpp"
#include "Stuff_Class.hpp"

void calc_non_deformation_forces(
    Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    Stuff_Class &stuff);