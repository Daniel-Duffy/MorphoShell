/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Calculate time step based on toy model stretching and bending analyses (take
whichever gives shortest characteristic time), and print.

*/


#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>


#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"
#include "calc_timestep.hpp"

void calc_timestep(
    Stuff_Class &stuff,
    Out_Stream_Class &log){

    const double pi = 3.14159265358979323846;

    // The 0.185 here came from some trial-and-improvement experimenting.
    // It's there so that if you redo a simulation with a coarser
    // mesh or smaller thickness such that the timestep switches from the bend
    // timestep to the stretch timestep, you shouldn't need to change timestep_prefactor
    // much (or hopefully at all) in the settings file.
    double stretch_timestep = 0.185 * stuff.timestep_prefactor * sqrt(stuff.ref_density / stuff.min_ref_shear_modulus) * stuff.approx_min_tri_size;
    double bend_timestep = stuff.timestep_prefactor * sqrt(stuff.ref_density*stuff.min_ref_thickness/stuff.char_flexural_rigidity)*stuff.approx_min_tri_size*stuff.approx_min_tri_size / (2*pi);

    log << "Short stretching and bending timescales: " << stretch_timestep << ", " << bend_timestep << std::endl;


    if( stuff.using_gradient_descent_dynamics ){
        stretch_timestep = 0.185 * stuff.timestep_prefactor * (stuff.damping_factor / (stuff.min_ref_shear_modulus*stuff.min_ref_thickness*stuff.min_ref_thickness)) * stuff.approx_min_tri_size * stuff.approx_min_tri_size / (2.0*pi);
        bend_timestep = stuff.timestep_prefactor * (stuff.damping_factor / (stuff.char_flexural_rigidity*stuff.min_ref_thickness)) * stuff.approx_min_tri_size * stuff.approx_min_tri_size * stuff.approx_min_tri_size * stuff.approx_min_tri_size / (2.0*pi * 2.0*pi * 2.0*pi);
    }


    if(stretch_timestep < bend_timestep){
        stuff.timestep = stretch_timestep;
        log << "Time step = " << stuff.timestep << ", set based on stretching rather than bending." << std::endl;
    }
    else{
        stuff.timestep = bend_timestep;
        log << "Time step = " << stuff.timestep << ", set based on bending rather than stretching." << std::endl;
    }

}
