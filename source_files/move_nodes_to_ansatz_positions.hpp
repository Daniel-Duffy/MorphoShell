//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "Stuff_Class.hpp"

void move_nodes_to_ansatz_positions(
    Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector< Eigen::Vector3d > &ansatz_node_positions,
    const Stuff_Class &stuff);