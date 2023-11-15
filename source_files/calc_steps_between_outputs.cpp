/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

Calculate steps_between_outputs based on dial_time / timestep, tuned by a
dimensionless parameter outputs_per_char_long_time in the setting file.
We do some rounding so steps_between_outputs is an int, but we ensure it's never 
zero; setting outputs_per_char_long_time to a huge number will give steps_between_outputs=1 instead, 
i.e. an output after every time step, which can be a useful thing sometimes. 
If outputs_per_char_long_time is set to a negative value we instead set steps_between_outputs 
to -1, which means no output files are regularly written at all.
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
#include "calc_steps_between_outputs.hpp"

void calc_steps_between_outputs(
    Stuff_Class &stuff,
    Out_Stream_Class &log){

    if( stuff.outputs_per_char_long_time < 0.0 ){
        stuff.steps_between_outputs = -1;
        stuff.output_regularly = false;
        log << "The setting printouts_per_char_long_time < 0, so regular printouts have been turned off completely." << std::endl;
    }
    else{
        if( stuff.outputs_per_char_long_time < 1.0 ){
            throw std::runtime_error("Error: To avoid potential divide-by-zero problems we do NOT allow stuff.printouts_per_char_long_time to lie in [0.0, 1.0). This is overkill somewhat, and could be relaxed if need be.");
        }
        else{
            stuff.output_regularly = true;
            stuff.steps_between_outputs = lround( stuff.char_long_time / (stuff.timestep * stuff.outputs_per_char_long_time) ) + 1;
            //log << "steps_between_outputs = " << stuff.steps_between_outputs << std::endl;
        }
    }


    // I usually set num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics
    // to either 0 or 10000, so the following lets me see 10 outputs from the gradient descent portion.
    if( stuff.using_gradient_descent_dynamics ){
        stuff.steps_between_outputs = 1000;
    }
}

