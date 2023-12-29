/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Function to advance the dynamics of the system by evolving node
positions and velocities according to the total forces on them. The scheme used
is so-called 'semi-implicit Euler' also called Symplectic Euler. This is first 
order, but it is simple and we are not interested in resolving the details of 
the dynamics, only the final equilibrium state. 

The forces are also checked for 'blowing up', which we treat as a 
kind of crash.
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

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <cmath>


#include "advance_dynamics.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void advance_dynamics(
    Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    Stuff_Class &stuff,
    Out_Stream_Class &log){
    
    // Check forces are well behaved. If not, throw error, and
    // display some stuff. Need to have this outside openMP #pragma 
    // because of openMP's exception handling rules.
    for(int d = 0; d < stuff.num_dofs; ++d){
        if( !(fabs(forces(d)) < 50000.0 * stuff.char_force_scale) ){
            log << " -----------------------------------------" << std::endl;
            log << " ------------CRASH REPORT-----------------" << std::endl;
            log << " -----------------------------------------" << std::endl;
            //log << "Offending node and its incident triangles: " << std::endl;
            log << "Offending dof is " << d << std::endl;
            log << "Offending force is " << forces(d) << std::endl;
            log << "Offending node is " << d%stuff.num_nodes << std::endl;
            //nodes[n].display();
            //for(int t = 0; t < nodes[n].incident_tri_ids.size(); ++t){
            //    triangles[nodes[n].incident_tri_ids(t)].display();
            //}
            throw std::runtime_error("implausibly_high_force");
        }
    }


    #pragma omp parallel for simd
    for(int d = 0; d < stuff.num_dofs; ++d){

        // Advance velocity before position (Semi-Implicit Euler, also
        // called Symplectic Euler).
        velocities(d) += (stuff.timestep / dof_masses(d)) * forces(d);

        // Gradient Descent dynamics may be used instead, overriding
        // the above velocity update. The damping and density bits here
        // come from thinking of gradient descent as the overdamped limit
        // of Newtonian dynamics.
        if( stuff.using_gradient_descent_dynamics ){
            velocities(d) = stuff.ref_density*stuff.min_ref_thickness*stuff.min_ref_thickness * forces(d) / (stuff.damping_factor * dof_masses(d));
        }
        
        // Advance position.
        dofs(d) += stuff.timestep * velocities(d);
    }


    // Advance upper 'glass slide'.
    stuff.upper_slide_z_coord += -stuff.timestep * stuff.slide_speed_prefactor * stuff.sample_char_length / stuff.bend_long_time;

    // Dial gravity, for a Calladine-esque setup. The minus sign in front of the speed prefactor is a hack so you can use
    // a single (negative) setting to have gravity increasing in strength while the top slide moves away without doing anything.
    // Using the slide_speed_prefactor as a gravitational dial speed prefactor is just a hack; adding a new setting would be better.
    //stuff.grav_field_strength += (-stuff.slide_speed_prefactor) * stuff.min_ref_shear_modulus * stuff.min_ref_thickness / (stuff.ref_density * stuff.total_ref_area);
}
