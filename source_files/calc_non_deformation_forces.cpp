/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

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
#include <Eigen/Sparse>
#include <omp.h>


#include "calc_non_deformation_forces.hpp"
#include "Node_Class.hpp"
#include "Stuff_Class.hpp"

void calc_non_deformation_forces(
    Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    Stuff_Class &stuff){

    #pragma omp parallel for simd
    for( int d = 0; d < stuff.num_dofs; ++d ){

        /* Linear viscous damping force to ensure equilibrium is approached.
        This damping model is chosen partly for simplicity but mainly for the
        desirable dispersion relation resulting from it.
        The factor of mass is chosen to counteract the damping overly penalising
        the boundary nodes, which will in general have lower masses.
        The damping force is only included for inertial Newtonian dynamics; it
        is not included in the 'force' vector for the alternative Gradient
        Descent dynamics. */
        if( !stuff.using_gradient_descent_dynamics ){
            forces(d) += -stuff.damping_factor * dof_masses(d) * velocities(d) / (stuff.ref_density*stuff.min_ref_thickness*stuff.min_ref_thickness);
        }


        // Gravity.
        //Eigen::Vector3d grav_field_direction(0.0, 0.0, -1.0);
        //if( d >= 2*stuff.num_nodes ){
        //    forces(d) += -stuff.grav_field_strength * dof_masses(d);
        //}
    }

    for( int n = 0; n < stuff.num_nodes; ++n ){

        // Squashing between two "glass slides". 
        // The lower slide is stationary, while 
        // the upper slide moves vertically (see 
        // advance_dynamics.cpp).
        double node_z_coord = dofs(2*stuff.num_nodes + n);
        double slide_z_coord;
        bool node_is_contacting_upper_slide = false;
        bool node_is_contacting_lower_slide = false;
        if( node_z_coord > stuff.upper_slide_z_coord ){
            slide_z_coord = stuff.upper_slide_z_coord;
            node_is_contacting_upper_slide = true;
        }
        if( node_z_coord < stuff.lower_slide_z_coord ){
            slide_z_coord = stuff.lower_slide_z_coord;
            node_is_contacting_lower_slide = true;
        }
        if( (node_is_contacting_upper_slide || node_is_contacting_lower_slide) && stuff.slide_stiffness_prefactor > 0.0 ){
            double slide_vert_force = stuff.slide_stiffness_prefactor * 
                                    stuff.min_ref_shear_modulus * 
                                    stuff.min_ref_thickness * 
                                    (slide_z_coord - node_z_coord);
            forces(2*stuff.num_nodes + n)  += slide_vert_force;
            if( node_is_contacting_upper_slide ){
                stuff.upper_slide_tot_vert_force += slide_vert_force;
            }
            if( node_is_contacting_lower_slide ){
                stuff.lower_slide_tot_vert_force += slide_vert_force;
            }

            // Slide friction. Sticking or sliding 
            // depending on the normal force. 
            // Simplest friction model I could think of.
            double non_friction_horiz_force_size = sqrt(forces(n)*forces(n) + forces(stuff.num_nodes + n)*forces(stuff.num_nodes + n));
            if( non_friction_horiz_force_size < stuff.slide_friction_coefficent * fabs(slide_vert_force) ){
                forces(n) = 0;
                forces(stuff.num_nodes + n) = 0;
                velocities(n) = 0;
                velocities(stuff.num_nodes + n) = 0;
            }
            else{
                forces(n) += -stuff.slide_friction_coefficent * fabs(slide_vert_force) * forces(n) / non_friction_horiz_force_size;
                forces(stuff.num_nodes + n) += -stuff.slide_friction_coefficent * fabs(slide_vert_force) * forces(stuff.num_nodes + n) / non_friction_horiz_force_size;
            }
        }
    }
}
