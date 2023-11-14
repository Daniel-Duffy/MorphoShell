
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

void calc_angle_deficits(
    Eigen::Matrix<double,Eigen::Dynamic,1> &angle_deficits,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    const std::vector<Node_Class> &nodes,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff);