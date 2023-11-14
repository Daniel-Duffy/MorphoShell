/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>


#include "Node_Class.hpp"
#include "Stuff_Class.hpp"

void perturb_node_positions_with_noise(
    Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions, 
    const Stuff_Class &stuff);