/*
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

Header file for set_up_patch_fitting.cpp function
*/

//Turn Eigen bounds checking off.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void set_up_patch_fitting(
    const std::vector<Node_Class> &nodes,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    Out_Stream_Class &log);
