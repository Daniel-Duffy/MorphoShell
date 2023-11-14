/* 
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this 
code in work that you publish or distribute.


/////////////////////////////////////////////////////

This file defines the member functions for the Node_Class class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there.*/

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
