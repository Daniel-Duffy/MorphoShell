
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "kahan_sum.hpp"

void calc_curvatures(
    Eigen::Matrix<double,Eigen::Dynamic,2> &curvatures,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &b_comps,
    const Stuff_Class &stuff);