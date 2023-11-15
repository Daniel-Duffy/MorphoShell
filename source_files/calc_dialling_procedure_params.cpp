/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

The dial factor will evolve from its starting value to 1.
The evolution is linear in time, except the dialling pauses
at a sequence of dial factor values, whose spacing is given 
by the setting step_between_dial_factors_to_hold_at. During these
pauses the simulation is in a 'holding'/'waiting' phase that allows the 
shell to come to equilibrium; once equilibrium is reached (to 
within the chosen tolerance) the code continues on to a new 
dialling phase. Each dialling phase takes dial_phase_time 
to complete, which equals the (approx) longest timescale 
in the simulation multiplied by the dial_phase_time_prefactor setting.
In each 'holding'/'waiting' phase, the check for equilibrium 
occurs at a rate determined by the time_between_equil_checks_prefactor setting.

When the dial_factor reaches 1, the final holding phase begins. Once it finishes
(equil is reached), the simulation terminates.

If the step_between_dial_factors_to_hold_at setting is > 1, then there will be 
only a single holding phase, which begins when the dial_factor reaches 1. The 
simulation terminates after that holding phase completes, as usual.
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
#include "calc_dialling_procedure_params.hpp"

void calc_dialling_procedure_params(
    Stuff_Class &stuff,
    Out_Stream_Class &log){

    if( (stuff.dial_factors_to_hold_at_spacing <= 0) || (stuff.dial_phase_time < 0) ){
        throw std::runtime_error(
            "\nstuff.step_between_dial_factors_to_hold_at and stuff.dial_phase_time must "
            "be >0 and >=0 respectively, and they aren't currently. Aborting.");
    }

    // Set stuff.dial_phase_time to be the (approx) longest timescale in the 
    // simulation multiplied by a prefactor chosen in the settings file.
    stuff.dial_phase_time = stuff.dial_phase_time_prefactor * stuff.char_long_time;
    stuff.timesteps_per_dial_phase = lround(stuff.dial_phase_time / stuff.timestep);
    log << "Duration of each dialling phase = " << stuff.dial_phase_time << " which is " << stuff.timesteps_per_dial_phase << " timesteps." << std::endl;

    // Calc how much the dial_factor needs to change per timestep (when dialling rather than waiting for equil).
    //stuff.dial_factor_change_per_timestep = stuff.step_between_dial_factors_to_hold_at * stuff.timestep / stuff.dial_phase_time;

    // Set stuff.timesteps_between_equil_checks using stuff.dial_phase_time multiplied 
    // by a prefactor chosen in the settings file.
    stuff.timesteps_between_equil_checks = lround(stuff.time_between_equil_checks_prefactor * stuff.dial_phase_time / stuff.timestep);
}


