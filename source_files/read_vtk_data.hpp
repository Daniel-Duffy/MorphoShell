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


#include <cstddef>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void read_vtk_data(
    std::vector<Node_Class> &nodes,
    Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    std::vector< Eigen::Vector3d > &ansatz_node_positions,
    std::vector<Tri_Class> &triangles,
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &abar_info,
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &bbar_info,
    Stuff_Class &stuff,
    double &dial_factor_to_start_from,
    Out_Stream_Class &log);
