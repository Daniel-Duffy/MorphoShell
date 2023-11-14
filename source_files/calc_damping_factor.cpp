/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Calculate 'dialling in' time and damping coefficient based on toy model
stretching and bending analyses. The 'Long times' are approximate characteristic
times for the longest-wavelength modes in the system for stretching and bending.
The damping coefficient is chosen to approximately critically damp the longest-
wavelength bending mode.
*/


#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>


#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"
#include "calc_damping_factor.hpp"

void calc_damping_factor(
    const double &damping_prefactor,
    Stuff_Class &stuff,
    Out_Stream_Class &log){

    const double pi = 3.14159265358979323846;

    const double minWavevector = 2.0 * pi / stuff.sample_char_length;
    const double damping_scale = 2.0 * sqrt(stuff.char_flexural_rigidity*stuff.ref_density*stuff.min_ref_thickness*stuff.min_ref_thickness*stuff.min_ref_thickness) * minWavevector * minWavevector;
    
    // Set the damping factor to be the characteristic damping scale for
    // critically damping long-wavelength modes multiplied by a supplied 
    // prefactor.
    stuff.damping_factor = damping_prefactor * damping_scale;
    log << "Numerical Damping Coefficient = " << stuff.damping_factor << std::endl;

}


