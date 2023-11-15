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
#include <Eigen/Sparse>

#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"


void calc_deformation_forces(
    Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    Eigen::Matrix<double,Eigen::Dynamic,1> &del_energy_by_del_continuum_quantities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &def_shear_moduli,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &def_thicknesses,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &bbar_comps,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_for_continuum_quantities_from_dofs_transpose,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    [[maybe_unused]] Out_Stream_Class &log);