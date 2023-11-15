/* 
This file defines the member functions for the Edge_Class class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there*/

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

#include <iostream>

#include "Edge_Class.hpp"


//This is a debugging tool to display the edge's data
void Edge_Class::display() {
    edge_log_stream.open();
    edge_log_stream << "-----------------------------" << std::boolalpha << std::endl;
    edge_log_stream << "Edge " << id << ":" << std::endl;
    //edge_log_stream << "Node ids " << node_ids.transpose() << std::endl;
    //edge_log_stream << "Adjacent triangle ids: " << adjacent_tri_ids.transpose() << std::endl;
    //edge_log_stream << "Left and Right 'other' (non-edge) node ids: " << otherNodeid_L << ", " << otherNodeid_R << std::endl;
    edge_log_stream << "Boundary indicator: " << is_boundary << std::endl;
    edge_log_stream << "-----------------------------" << std::endl;
    edge_log_stream.close();
}
