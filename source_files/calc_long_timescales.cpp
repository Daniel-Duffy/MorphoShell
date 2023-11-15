/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Calculate 'dialling in' time and damping coefficient based on toy model
stretching and bending analyses. The 'Long times' are approximate characteristic
times for the longest-wavelength modes in the system for stretching and bending.
The damping coefficient is chosen to approximately critically damp the longest-
wavelength bending mode.
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
#include <vector>
#include <cmath>
#include <stdexcept>


#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"
#include "calc_long_timescales.hpp"

void calc_long_timescales(
    Stuff_Class &stuff,
    Out_Stream_Class &log){

    const double pi = 3.14159265358979323846;

    const double approx_min_wavevector = 2.0 * pi / stuff.sample_char_length;

    stuff.stretch_long_time = 2.0 * pi * sqrt(stuff.ref_density / stuff.min_ref_shear_modulus) / approx_min_wavevector;
    stuff.bend_long_time =    2.0 * pi * sqrt(stuff.ref_density * stuff.min_ref_thickness / stuff.char_flexural_rigidity) / (approx_min_wavevector * approx_min_wavevector);

    if(stuff.stretch_long_time > stuff.bend_long_time){
        stuff.char_long_time = stuff.stretch_long_time;    
    }
    else{
        stuff.char_long_time = stuff.bend_long_time;
    }

    log << "The long characteristic time = " << stuff.char_long_time << std::endl;
}


