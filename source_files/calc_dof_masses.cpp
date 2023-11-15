/* 
Function to finish setting up initial data for the flat LCE sheet: setting
node velocities to zero, and storing the initial in-plane sides' components for
the triangles (using the x-y plane basis), which will then not change. The
initial areas are also calculated and stored. Also, calculate node masses by
having each triangle contribute 1/3 of its initial mass to each of its vertcies.
*/

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
#include <cmath>

#include "calc_dof_masses.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void calc_dof_masses(
    Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    [[maybe_unused]] Out_Stream_Class &log){


    for(int d = 0; d < stuff.num_dofs; ++d){
        //Set all nodes masses to zero before calculating them next
        dof_masses(d) = 0.0;
    }

    for( const auto& tri: triangles ){

        /* Add 1/3 of the mass of each triangle to each of its vertices.*/
        for(int v = 0; v < 3; ++v){
            dof_masses[tri.vertex_ids[v]] += stuff.ref_density * tri.ref_area * tri.ref_thickness / 3.0;
            dof_masses[tri.vertex_ids[v]+stuff.num_nodes] += stuff.ref_density * tri.ref_area * tri.ref_thickness / 3.0;
            dof_masses[tri.vertex_ids[v]+2*stuff.num_nodes] += stuff.ref_density * tri.ref_area * tri.ref_thickness / 3.0;
        }
    }
}




















