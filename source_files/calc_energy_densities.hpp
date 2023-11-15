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

#include "Stuff_Class.hpp"

void calc_energy_densities(
    Eigen::Matrix<double,Eigen::Dynamic,2> &energy_densities,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &a_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &b_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &bbar_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &def_shear_moduli,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &def_thicknesses,
    const Stuff_Class &stuff,
    [[maybe_unused]] Out_Stream_Class &log);