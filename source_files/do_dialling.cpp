/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/

//Turn Eigen bounds checking off for speed.
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
    const Eigen::Matrix<double,Eigen::Dynamic,3> &initial_a_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &initial_b_comps,
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
    // happens in two parts: one here, and one at the end of this
    // file.
    double dial_factor_copy = dial_factor;
    if( stuff.dialling_from_ansatz_rather_than_ref_state ){
        // This is so the subsequent dialling protocols set
        // abar and bbar to their final values.
        dial_factor = 1.0;
    }


    // Uncomment the desired dialling protocol below (remember
    // the abar_info and bbar_info you provided in the input file
    // must hold the correct info for your desired protocol!). 
    // Feel free to add your own!
    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////

    // This case assumes that abar_info and bbar_info 
    // are just equal to the final (dial_factor=1) 
    // abar_comps and bbar_comps, and we just dial 
    // abar and bbar linearly from their initial
    // values to those final values.
    
    abar_comps = (1.0-dial_factor) * initial_a_comps + dial_factor * abar_info;
    bbar_comps = (1.0-dial_factor) * initial_b_comps + dial_factor * bbar_info;
    
    /////////////////////////////////////////////////
    /*
    // As just above, except we dial bbar first completely
    // (reaching bbar_info at dial_factor = 0.5) and then 
    // dial abar completely.
    if( dial_factor < 0.5 ){
        abar_comps = initial_a_comps;
        bbar_comps = (1.0-2*dial_factor) * initial_b_comps + 2*dial_factor * bbar_info;
    }
    else{
        abar_comps = (1.0-2*(dial_factor-0.5)) * initial_a_comps + 2*(dial_factor-0.5) * abar_info;
        bbar_comps = bbar_info;
    }
    */

    /////////////////////////////////////////////////

    // This case assumes that bbar = 0 and that
    // abar is of LCE form, and we dial the LCE
    // lambda linearly in time from 0. So this 
    // is ideal for LCEs with in-plane director 
    // patterns. Each row of abar_info should hold 
    // director angle, final (target) lambda value, opto-thermal poisson ratio.
    /*
    bbar_comps.fill(0.0);
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        double angle = abar_info(t,0);
        double lambda = (1.0-dial_factor) * 1.0 + dial_factor * abar_info(t,1);
        double nu_ot = abar_info(t,2);

        abar_comps(t,0) = lambda * lambda * cos(angle) * cos(angle)  +  pow(lambda, -2.0 * nu_ot) * sin(angle) * sin(angle);
        abar_comps(t,1) = (lambda * lambda - pow(lambda, -2.0 * nu_ot)) * cos(angle) * sin(angle);
        abar_comps(t,2) = pow(lambda, -2.0 * nu_ot) * cos(angle) * cos(angle)  +  lambda * lambda * sin(angle) * sin(angle);
    }
    */

    /////////////////////////////////////////////////

    // This case assumes that bbar = 0 and that
    // abar is of LCE form, and we dial the LCE
    // lambda linearly in time from one inhomogeneous
    // lambda field to the next --- there are 11,
    // corresponding to dial_factor = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 .
    // We must therefore have dial_factors_to_hold_at_spacing = 0.1 in settings.
    /*
    bbar_comps.fill(0.0);
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){
                
        double nu_ot = abar_info(t,11+1-1);
        double angle = abar_info(t,11+1);
        int prev_dial_factor_checkpoint_idx = static_cast<int>(std::floor(dial_factor/stuff.dial_factors_to_hold_at_spacing)); // Fine as long as argument of cast is not large, in which case might have a double outside the range of possible ints.
        double fraction_of_way_through_current_dialling_phase = (dial_factor/ stuff.dial_factors_to_hold_at_spacing-prev_dial_factor_checkpoint_idx);
        double prev_dial_factor_checkpoint = abar_info(t, prev_dial_factor_checkpoint_idx);
        double next_dial_factor_checkpoint = abar_info(t, prev_dial_factor_checkpoint_idx+1);
        double lambda = (1-fraction_of_way_through_current_dialling_phase) * prev_dial_factor_checkpoint + fraction_of_way_through_current_dialling_phase * next_dial_factor_checkpoint;

        //if( t == 0 ){log << "lambda  " << lambda << "  frac  " << fraction_of_way_through_current_dialling_phase << "   idx   " << prev_dial_factor_checkpoint_idx << std::endl;}
        

        abar_comps(t,0) = lambda * lambda * cos(angle) * cos(angle)  +  pow(lambda, -2.0 * nu_ot) * sin(angle) * sin(angle);
        abar_comps(t,1) = (lambda * lambda - pow(lambda, -2.0 * nu_ot)) * cos(angle) * sin(angle);
        abar_comps(t,2) = pow(lambda, -2.0 * nu_ot) * cos(angle) * cos(angle)  +  lambda * lambda * sin(angle) * sin(angle);
    }
    */

    ////////////////////////////////////////////////

    // You could do an example here of dialling and LCE 
    // except with lambda modified by height of the tri
    // centroid or angle of normal to the z axis.


    /////////////////////////////////////////////////

    // This case assumes that bbar = 0 and uses isotropic
    // swelling to encode constant negative
    // Gauss curvature K, with fabs(K) = 2 * dial_factor.
    // The r values in abar_info are just values of radial
    // coordinate in the 2D reference plane.
    /*
    bbar_comps.fill(0.0);
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        double r = abar_info(t,0);

        double absK = dial_factor * 2.0;

        double Omega = 1.0 / (1.0 - 0.25 * r * r * absK); // Lineal swelling factor

        abar_comps(t,0) = Omega * Omega;
        abar_comps(t,1) = 0;
        abar_comps(t,2) = Omega * Omega;
    }
    */


    ////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    // Now we do the second part of the dialling-from-ansatz
    // procedure that we began at the start of this file.
    // In general watch out for "aliasing" with Eigen matrices,
    // though here we are fine.
    if( stuff.dialling_from_ansatz_rather_than_ref_state ){

        // Usual dialling from ansatz.
        
        abar_comps = (1.0-dial_factor_copy) * initial_a_comps + dial_factor_copy * abar_comps;
        bbar_comps = (1.0-dial_factor_copy) * initial_b_comps + dial_factor_copy * bbar_comps;
        

        // Or instead do no dialling so the ansatz abar and bbar are just abar and bbar for all time.
        /*
        abar_comps = initial_a_comps;
        bbar_comps = initial_b_comps;
        */

        // Or for Mingchao's flying saucer, abar and bbar just stay put except in a central reference disk where they
        // decrease, corresponding to a shrinkage by Omega and a programmed-curvature increase of 1/Omega (the
        // shrinking piece was a sphere.)
        /*
        abar_comps = initial_a_comps;
        bbar_comps = initial_b_comps;
        // Nodes corresponding to shrinking region will be tagged in python.
        for(int t = 0; t < stuff.num_tris; ++t){
            if( triangles[t].tag ){
                abar_comps.row(t) = (1.0-dial_factor_copy)*(1.0-dial_factor_copy) * initial_a_comps.row(t);
                bbar_comps.row(t) = (1.0-dial_factor_copy) * initial_b_comps.row(t);
            }
        }
        */
    }




    // NOW LET'S DIAL DEF_THICKNESS.
    // UNCOMMENT THE APPROPRIATE LINE, OR 
    // ADD YOUR OWN!
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // No changing of thickness during deformation.
        def_thicknesses(t) = triangles[t].ref_thickness; 

        // LCE with "opto-thermal Poisson ratio" of 1/2.
        //def_thicknesses(t) = triangles[t].ref_thickness / sqrt(abar_comps(t,0)*abar_comps(t,2) - abar_comps(t,1)*abar_comps(t,1));

    }

    // NOW LET'S DIAL DEF_SHEAR_MODULUS.
    // UNCOMMENT THE APPROPRIATE LINE, OR 
    // ADD YOUR OWN!
    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // No changing of shear modulus during deformation.
        def_shear_moduli(t) = triangles[t].ref_shear_modulus; 

    }
}
