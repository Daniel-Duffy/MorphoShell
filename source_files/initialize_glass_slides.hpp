
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <string>

#include "Stuff_Class.hpp"
#include "find_extremum_of_doubles.hpp"


void initialize_glass_slides(
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    Stuff_Class &stuff);