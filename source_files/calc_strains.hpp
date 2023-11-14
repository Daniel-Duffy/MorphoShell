
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "Stuff_Class.hpp"

void calc_strains(
    Eigen::Matrix<double,Eigen::Dynamic,1> &strains,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &a_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    const Stuff_Class &stuff);