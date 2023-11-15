/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to advance the dynamics of the system by evolving node
positions and velocities according to the total forces on them. The scheme used
is so-called 'semi-implicit Euler' also called Symplectic Euler. This is first 
order, but it is simple and we are not interested in resolving the details of 
the dynamics, only the final equilibrium state. 

The forces are also checked, to catch code crashed in which the forces usually
'blow up'.
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
#include <iostream>
#include <cmath>
#include <omp.h>


#include "impose_constraints.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void impose_constraints(
    const std::vector<Node_Class> &nodes,
    [[maybe_unused]] Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    [[maybe_unused]] Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    [[maybe_unused]] Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    [[maybe_unused]] const std::vector<Tri_Class> &triangles,
    [[maybe_unused]] const Stuff_Class &stuff){

    #pragma omp parallel for simd
    for( int n = 0; n < stuff.num_nodes; ++n ){

        if( nodes[n].is_constrained ){

            // Full clamp. Can comment
            // out lines to e.g. clamp
            // only in z.
            forces(n) = 0.0;
            forces(n+stuff.num_nodes) = 0.0;
            forces(n+2*stuff.num_nodes) = 0.0;
        }
    }
}
