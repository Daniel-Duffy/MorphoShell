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

*/


//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>

#include <cstddef>
#include <string>
#include <fstream>
#include <iomanip> //for setting output precision etc
#include <vector>
#include <stdexcept>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"

void write_output(
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    const std::vector<Tri_Class> &triangles,
    const int &stepcount,
    const double &dial_factor,
    const Eigen::Matrix<double,Eigen::Dynamic,2> &curvatures,
    const Eigen::Matrix<double,Eigen::Dynamic,2> &energy_densities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &angle_deficits,
    const Eigen::Matrix<double,Eigen::Dynamic,8> &stresses,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &strains,
    [[maybe_unused]] const std::vector<double> &energies,
    const Stuff_Class &stuff);
