/* 
This file defines the member functions for the Node_Class class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there.*/

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
#include <iomanip>

#include "Node_Class.hpp"

//This is a debugging tool to display the node's data
/*
void Node_Class::display(){
    node_log_stream.open();
    node_log_stream << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    node_log_stream << "Node " << id << ":" << std::endl;
    //node_log_stream << "ids of incident triangles: " << "\n" << incident_tri_ids << std::endl;
    //node_log_stream << "neighbourNodeids = " << "\n" << neighbour_node_ids << std::endl;
    node_log_stream << "Position = " << "\n" << pos << std::endl;
    node_log_stream << "Velocity = " << "\n" << vel << std::endl;
    node_log_stream << "Force = " << "\n" << force << std::endl;
    node_log_stream << "Mass = " << mass << std::endl;
    node_log_stream << "Boundary indicator: " << is_boundary << std::endl;
    node_log_stream << "Clamp indicator: " << clamp_indicator << std::endl;
    node_log_stream << "Load indicator: " << load_indicator << std::endl;
    node_log_stream << "-----------------------------" << std::endl;
    node_log_stream.close();
}
*/
