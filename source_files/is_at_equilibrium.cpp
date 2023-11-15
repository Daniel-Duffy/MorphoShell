/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to determine whether the shell is in equilibrium or not (to within
chosen thresholds), based on:
1) The ratio of the largest node force to a typical simulation force scale given
by shear modulus * smallest initial element size * SheetThickness.
2) The largest value of the ratio of node velocity to a typical local velocity
scale - the (stretching) sound speed.

Returns true if equilibrium has been reached (to within the chosen tolerance),
and false otherwise.

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

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <omp.h>


#include "is_at_equilibrium.hpp"
#include "Node_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

#include "find_extremum_of_doubles.hpp"


bool is_at_equilibrium(
    const Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Stuff_Class &stuff,
    Out_Stream_Class &log){

    // These scalar forces are norms of force vectors.
    std::vector<double> force_norms(stuff.num_nodes);
    std::vector<double> node_speeds(stuff.num_nodes);

    // Calculate node non-damping forces and scaled speeds.
    #pragma omp parallel for simd
    for( int n = 0; n < stuff.num_nodes; ++n ){

        Eigen::Vector3d force_vec = {forces(n), forces(n+stuff.num_nodes), forces(n+2*stuff.num_nodes)};
        Eigen::Vector3d vel_vec = {velocities(n), velocities(n+stuff.num_nodes), velocities(n+2*stuff.num_nodes)};

        // In a way the nicest thing here is to look at norms of just the non-damping part of the force, since that's the but
        // that tells you how far from the minimum of elastic potential energy you are. However it's a bit fiddly to do because
        // we've already combined damping and non-damping forces together. The condition that node speeds are small puts a bound
        // on the damping forces though, so as long as both thresholds are sufficiently small, you still end up being sure that
        // non-damping forces and speeds are both small, which is all that really matters.
        force_norms[n] = force_vec.norm();
        //non_damping_forces[n] = (force_vec + (stuff.damping_factor * dof_masses(n) * vel_vec / stuff.ref_density)).norm();
        node_speeds[n] = vel_vec.norm();
    }

    double max_force = min_or_max_of_doubles(force_norms, std::string("max"));
    size_t max_force_node = argmin_or_argmax_of_doubles(force_norms, std::string("max"));

    // By scaled speed we mean node speed divided 
    // by local stretching sound (wave) speed.
    double max_speed = min_or_max_of_doubles(node_speeds, std::string("max"));
    size_t max_speed_node = argmin_or_argmax_of_doubles(node_speeds, std::string("max"));


    /* Now check whether non-damping force and speed ratios are below chosen
    `equilibrium' thresholds.*/
    log << "Max ratio of non-damping force to characteristic force = " << max_force/stuff.char_force_scale << " (node " << max_force_node << ")" << "\n";
    log << "Max ratio of node speed to characteristic speed = " << max_speed/stuff.char_speed_scale << " (node " << max_speed_node << ")" << "\n";

    if( (max_force / stuff.char_force_scale) < stuff.force_to_char_force_ratio_equil_threshold
     && (max_speed/stuff.char_speed_scale) < stuff.speed_to_char_speed_ratio_equil_threshold ){
        return true;
    }
    else{
        return false;
    }
}
