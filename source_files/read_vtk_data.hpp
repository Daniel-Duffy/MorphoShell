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

Header file for readVTKData.cpp function.
*/

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
