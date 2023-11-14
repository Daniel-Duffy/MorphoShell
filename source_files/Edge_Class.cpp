/* 
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this 
code in work that you publish or distribute.

This file is part of Shellmorph.

Shellmorph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Shellmorph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

This file defines the member functions for the Edge_Class class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there*/

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
