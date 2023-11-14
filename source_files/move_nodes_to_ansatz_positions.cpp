/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "move_nodes_to_ansatz_positions.hpp"
#include "Stuff_Class.hpp"

void move_nodes_to_ansatz_positions(
    Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector< Eigen::Vector3d > &ansatz_node_positions,
    const Stuff_Class &stuff){

    for(int n = 0; n < stuff.num_nodes; ++n){
        node_positions(n,0) = ansatz_node_positions.at(n)(0);
        node_positions(n,1) = ansatz_node_positions.at(n)(1);
        node_positions(n,2) = ansatz_node_positions.at(n)(2);
    }

    // Delete ansatz_node_positions (no longer needed).
    ansatz_node_positions.resize(0);
}
