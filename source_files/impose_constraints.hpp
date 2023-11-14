//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void impose_constraints(
    const std::vector<Node_Class> &nodes,
    [[maybe_unused]] Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    [[maybe_unused]] Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    [[maybe_unused]] Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    [[maybe_unused]] const std::vector<Tri_Class> &triangles,
    [[maybe_unused]] const Stuff_Class &stuff);