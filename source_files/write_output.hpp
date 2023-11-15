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
