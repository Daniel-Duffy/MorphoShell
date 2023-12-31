/////////////////////////////////////////////////////
/*
Copyright (C) 2023, Daniel Duffy, daniellouisduffy@gmail.com. All rights reserved.
Please cite Daniel Duffy and John S. Biggins if you 
use any part of this code in work that you publish or distribute.

This file is part of MorphoShell.

MorphoShell is distributed under the terms of the Cambridge Academic
Software License (CASL). You should have received a copy of the license
along with MorphoShell. If not, contact Daniel Duffy, daniellouisduffy@gmail.com.
*/
/////////////////////////////////////////////////////

// Turn Eigen bounds checking off for speed.
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
