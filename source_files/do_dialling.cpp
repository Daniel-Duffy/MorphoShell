/*
Here we 'dial' various quantities, most notably abar and bbar.
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

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <cstddef>

#include "do_dialling.hpp"
#include "Tri_Class.hpp"
#include "Node_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void do_dialling(
    Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    Eigen::Matrix<double,Eigen::Dynamic,3> &bbar_comps,
    Eigen::Matrix<double,Eigen::Dynamic,1> &def_thicknesses,
    Eigen::Matrix<double,Eigen::Dynamic,1> &def_shear_moduli,
    double dial_factor,
    const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &abar_info,
    [[maybe_unused]] const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &bbar_info,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &a_comps_at_zero_dial_factor,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &b_comps_at_zero_dial_factor,
    [[maybe_unused]] const std::vector<Tri_Class> &triangles,
    [[maybe_unused]] const std::vector<Node_Class> &nodes,
    [[maybe_unused]] const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    [[maybe_unused]] const Eigen::Matrix<double,Eigen::Dynamic,3> &normals,
    [[maybe_unused]] Out_Stream_Class &log,
    const Stuff_Class &stuff){


    // FIRST LET'S DIAL ABAR AND BBAR.

    // If you chose to dial abar and bbar from the ansatz
    // a and b values, then abar and bbar are dialled linearly
    // from the ansatz a and b to the final (dial_factor=1)
    // abar and bbar. This is a bit fiddly to implement, and
    // happens in two parts: one here, and later on in this
    // file.
    double dial_factor_copy = dial_factor;
    if( stuff.dialling_from_ansatz_rather_than_ref_state ){
        // This is so the subsequent dialling protocols set
        // abar and bbar to their final values.
        dial_factor = 1.0;
    }


    // Uncomment (only) the desired dialling protocol below (remember
    // the abar_info and bbar_info you provided in the input file
    // must hold the correct info for your desired protocol!). 
    // Feel free to add your own!
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    
    // 1)
    // This case assumes that abar_info and bbar_info 
    // are just equal to the final (dial_factor=1) 
    // abar_comps and bbar_comps, and we just dial 
    // abar and bbar linearly from their initial
    // values to those final values.
    //abar_comps = (1.0-dial_factor) * a_comps_at_zero_dial_factor + dial_factor * abar_info;
    //bbar_comps = (1.0-dial_factor) * b_comps_at_zero_dial_factor + dial_factor * bbar_info;
    
    
    /////////////////////////////////////////////////
    /*
    // 2)
    // As just above, except we dial bbar first completely
    // (reaching bbar_info at dial_factor = 0.5) and then 
    // dial abar completely.
    if( dial_factor < 0.5 ){
        abar_comps = a_comps_at_zero_dial_factor;
        bbar_comps = (1.0-2*dial_factor) * b_comps_at_zero_dial_factor + 2*dial_factor * bbar_info;
    }
    else{
        abar_comps = (1.0-2*(dial_factor-0.5)) * a_comps_at_zero_dial_factor + 2*(dial_factor-0.5) * abar_info;
        bbar_comps = bbar_info;
    }
    */

    /////////////////////////////////////////////////
    /*
    // 3)
    // This case assumes that abar is of LCE form for an in-plane director, 
    // (contraction by lambda along the director and lambda^(-nu_therm) orthogonally), 
    // and we dial the LCE lambda linearly in time from 1 to its final value. So this 
    // is ideal for LCEs with in-plane director 
    // patterns. Each row of abar_info should hold 
    // director angle, final (target) lambda value, opto-thermal poisson ratio.
    // Here we also assume bbar = 0.
    bbar_comps.fill(0.0);
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        double angle = abar_info(t,0); // Director angle relative to reference-state x-axis.
        double lambda = (1.0-dial_factor) * 1.0 + dial_factor * abar_info(t,1);
        double nu_ot = abar_info(t,2);

        abar_comps(t,0) = lambda * lambda * cos(angle) * cos(angle)  +  pow(lambda, -2.0 * nu_ot) * sin(angle) * sin(angle);
        abar_comps(t,1) = (lambda * lambda - pow(lambda, -2.0 * nu_ot)) * cos(angle) * sin(angle);
        abar_comps(t,2) = pow(lambda, -2.0 * nu_ot) * cos(angle) * cos(angle)  +  lambda * lambda * sin(angle) * sin(angle);
    }
    */

        /////////////////////////////////////////////////
    
    // 4)
    // This case assumes that abar is of LCE form for an in-plane director, 
    // (contraction by lambda along the director and lambda^(-nu_therm) orthogonally), 
    // and we dial the LCE lambda linearly in time from 1 to its final value. So this 
    // is ideal for LCEs with in-plane director 
    // patterns. Each row of abar_info should hold 
    // director angle, final (target) lambda value, opto-thermal poisson ratio.
    // Here we also assume that there is a preferred curvature (bbar), specifically with
    // the preferred principal curvature directions lying parallel to the director and orthogonal to it.
    // This might (for example) correspond to a setup in which the top of the sheet contracts by more 
    // than the bottom of the sheet (but along the same director on both sides). Thus there are two 
    // preferred principal curvatures, kappabar_par and kappabar_perp. Here both kappabar_par and 
    // kappabar_perp are dialled linearly over time, from 0 to their final target values (given in bbar_info).
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        double angle = abar_info(t,0); // Director angle relative to reference-state x-axis.
        double lambda = (1.0-dial_factor) * 1.0 + dial_factor * abar_info(t,1);
        double nu_ot = abar_info(t,2);

        abar_comps(t,0) = lambda * lambda * cos(angle) * cos(angle)  +  pow(lambda, -2.0 * nu_ot) * sin(angle) * sin(angle);
        abar_comps(t,1) = (lambda * lambda - pow(lambda, -2.0 * nu_ot)) * cos(angle) * sin(angle);
        abar_comps(t,2) = pow(lambda, -2.0 * nu_ot) * cos(angle) * cos(angle)  +  lambda * lambda * sin(angle) * sin(angle);
        
        //double angle = bbar_info(t,0); // Director angle relative to reference-state x-axis.
        double kappabar_par = dial_factor * bbar_info(t,1);
        double kappabar_perp = dial_factor * bbar_info(t,2);

        // Eigenvalues of matrix of downstairs components of bbar.
        double bb1 = lambda * lambda * kappabar_par;
        double bb2 = pow(lambda, -2.0 * nu_ot) * kappabar_perp;

        bbar_comps(t,0) = bb1 * cos(angle) * cos(angle)  +  bb2 * sin(angle) * sin(angle);
        bbar_comps(t,1) = (bb1 - bb2) * cos(angle) * sin(angle);
        bbar_comps(t,2) = bb2 * cos(angle) * cos(angle)  +  bb1 * sin(angle) * sin(angle);

        // Now the programmed shape operator Inverse[abar].bbar has eigenvalues kappabar_par and kappabar_perp,
        // and corresponding eigenvectors that are the director and a vector orthogonal to the director.
    }
    

    ////////////////////////////////////////////////
    // 5)
    // Could include an example here of dialling and LCE 
    // except with lambda modified by height of the tri
    // centroid or angle of normal to the z axis.



    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Now we do the second part of the dialling-from-ansatz
    // procedure that we began at the start of this file.
    // In general watch out for "aliasing" with Eigen matrices,
    // though here we are fine.
    if( stuff.dialling_from_ansatz_rather_than_ref_state ){

        // Usual dialling from ansatz.
        // A)
        abar_comps = (1.0-dial_factor_copy) * a_comps_at_zero_dial_factor + dial_factor_copy * abar_comps;
        bbar_comps = (1.0-dial_factor_copy) * b_comps_at_zero_dial_factor + dial_factor_copy * bbar_comps;
        
        // B)
        // Or instead do no dialling so the ansatz a and b are just abar and bbar for all time.
        /*
        abar_comps = a_comps_at_zero_dial_factor;
        bbar_comps = b_comps_at_zero_dial_factor;
        */

        // C)
        // Or for Mingchao's flying saucer, abar and bbar just stay put except in a central reference disk where they
        // decrease, corresponding to a shrinkage by Omega and a programmed-curvature increase of 1/Omega (the
        // shrinking piece was a sphere.)
        /*
        abar_comps = a_comps_at_zero_dial_factor;
        bbar_comps = b_comps_at_zero_dial_factor;
        // Nodes corresponding to shrinking region will be tagged in python.
        for(int t = 0; t < stuff.num_tris; ++t){
            if( triangles[t].tag ){
                abar_comps.row(t) = (1.0-dial_factor_copy)*(1.0-dial_factor_copy) * a_comps_at_zero_dial_factor.row(t);
                bbar_comps.row(t) = (1.0-dial_factor_copy) * b_comps_at_zero_dial_factor.row(t);
            }
        }
        */
        
        // D)
        // Keep abar as the ansatz a for all time, but dial bbar in the usual manner.
        //abar_comps = a_comps_at_zero_dial_factor;
        //bbar_comps = (1.0-dial_factor_copy) * b_comps_at_zero_dial_factor + dial_factor_copy * bbar_comps;
    }




    // NOW LET'S DIAL DEF_THICKNESS.
    // UNCOMMENT THE APPROPRIATE LINE, OR 
    // ADD YOUR OWN!
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // 1)
        // No changing of thickness during deformation.
        def_thicknesses(t) = triangles[t].ref_thickness; 

        // 2)
        // LCE with "opto-thermal Poisson ratio" of 1/2.
        //def_thicknesses(t) = triangles[t].ref_thickness / sqrt(abar_comps(t,0)*abar_comps(t,2) - abar_comps(t,1)*abar_comps(t,1));

    }

    // NOW LET'S DIAL DEF_SHEAR_MODULUS.
    // UNCOMMENT THE APPROPRIATE LINE, OR 
    // ADD YOUR OWN!
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // 1)
        // No changing of shear modulus during deformation.
        def_shear_moduli(t) = triangles[t].ref_shear_modulus; 

    }
}
