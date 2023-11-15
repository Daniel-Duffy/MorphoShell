// Main file for MorphoShell.

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

#include <iostream> // Used for outputting messages to the terminal
#include <assert.h> // Used for debugging tools
#include <stdexcept> // Used for exception (error) handling
#include <cmath> // Used for some simple maths functions
#include <Eigen/Dense> // Used for matrices of numbers
#include <Eigen/Sparse>
#include <string> // Used for creating directory for output data
#include <vector>// Used for some vectors
#include <iomanip> // For setting time and date format, and std::out precision if needed

#include "get_real_time.hpp"
#include "extract_filename.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"
#include "preliminary_setup.hpp"
#include "read_settings.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Edge_Class.hpp"
#include "read_vtk_data.hpp"
#include "calc_mesh_info.hpp"
#include "set_up_patch_fitting.hpp"
#include "calc_deriv_mats.hpp"
#include "calc_long_timescales.hpp"
#include "calc_damping_factor.hpp"
#include "calc_timestep.hpp"
#include "calc_steps_between_outputs.hpp"
#include "calc_char_mechanics_scales.hpp"
#include "Status_Enum.hpp"
#include "calc_dialling_procedure_params.hpp"
#include "perturb_node_positions_with_noise.hpp"
#include "calc_dof_masses.hpp"
#include "move_nodes_to_ansatz_positions.hpp"
#include "initialize_glass_slides.hpp"
#include "calc_surface_derivs.hpp"
#include "calc_a_comps_and_b_comps_and_normals.hpp"
#include "do_dialling.hpp"
#include "is_at_equilibrium.hpp"
#include "calc_energies.hpp"
#include "write_output.hpp"
#include "calc_curvatures.hpp"
#include "calc_angle_deficits.hpp"
#include "calc_stresses.hpp"
#include "calc_strains.hpp"
#include "impose_constraints.hpp"
#include "advance_dynamics.hpp"
#include "calc_deformation_forces.hpp"
#include "calc_non_deformation_forces.hpp"
#include "calc_energy_densities.hpp"



int main(int argc, char *argv[]){

// First, some checks and setting
// up the output directory, log file, etc.
Out_Stream_Class log;
Stuff_Class stuff; // This will hold settings and other useful quantities.
try{ 
preliminary_setup(
    argc, 
    argv,
    stuff,
    log);} catch(const std::runtime_error &error){std::cerr << error.what() << std::endl; return -1;}


// Now we read in the settings file.
// Settings are stored in 'stuff', which
// will also contain other useful quantities.
try{ read_settings(stuff, log); } 
catch(const std::out_of_range &out_of_bounds_error){
    log << "Vector out-of-bounds error while reading settings: " << out_of_bounds_error.what() <<
     "\nFirst thing to check: are all the settings read_settings.cpp is trying to read in your settings file?" << std::endl; return -1;}                             
catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}


// Now we read our input data.
std::vector<Node_Class> nodes;
Eigen::Matrix<double,Eigen::Dynamic,3> node_positions;
std::vector< Eigen::Vector3d > ansatz_node_positions;
std::vector<Tri_Class> triangles;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>  abar_info;
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> bbar_info;
double dial_factor_to_start_from = 0.0;
try{ read_vtk_data(
    nodes,
    node_positions,
    ansatz_node_positions, 
    triangles, 
    abar_info, 
    bbar_info, 
    stuff,  
    dial_factor_to_start_from, 
    log);} 
    catch(const std::out_of_range &out_of_bounds_error){log << "Vector out-of-bounds error while reading vtk: " << out_of_bounds_error.what() << std::endl; return -1;}
    catch(const std::runtime_error &error){
        log << error.what() << std::endl;
        log << "An error while reading data likely implies a problem with a data \n" 
        "file, for example a mismatch between the \n"
        "number of nodes or triangles stated and the number actually present; \n"
        "or other similar mismatches in numbers of data values; or a format problem. \n"
        "Remember the input file format is stricter than standard vtk." << std::endl;
        return -1;}
if( stuff.dialling_from_ansatz_rather_than_ref_state ){
    dial_factor_to_start_from = 0.0;
}



// Print some things.
log << "Dial factor will start from " << dial_factor_to_start_from << std::endl;
log << "Number of nodes = " << stuff.num_nodes << std::endl;
log << "Number of triangles = " << stuff.num_tris << std::endl;
// Print warning if number of triangles is low - partly because
// properly resolving out features of the problem is important for 
// physical accuracy, and partly to possible bugs in this rather special case.
if( stuff.num_tris < 50 ){
    log << "Your mesh has a small number of triangles. "
    "This is probably not a good idea for physics, "
    "and also bugs are more likely in extreme cases." << std::endl;
}


// Now we calculate some more information about our mesh.
std::vector<Edge_Class> edges;
try{ calc_mesh_info(
    nodes, 
    node_positions, 
    triangles, 
    edges, 
    stuff, 
    log);} catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}

// Create vector of dofs.
Eigen::Matrix<double,Eigen::Dynamic,1> dofs(stuff.num_dofs, 1);


// Set up patch-fitting approach used to approx surface second derivatives.
try {set_up_patch_fitting( 
    nodes, 
    node_positions, 
    triangles, 
    stuff, 
    log);} catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}


// Set up matrices that multiply onto dofs to give various derivatives of the surface.
// The reason we store the transpose separately (very wasteful in terms of memory) is that we
// Eigen will not multi-thread a ColMajor-sparse * dense multiplication, and the transpose of 
// a RowMajor-sparse is effectively ColMajor. Storing the transpose as RowMajor allows Eigen to
// multithread when we use it, which is good for speed.
Eigen::SparseMatrix<double, Eigen::RowMajor> mat_for_continuum_quantities_from_dofs(stuff.num_continuum_quantities, stuff.num_dofs);
Eigen::SparseMatrix<double, Eigen::RowMajor> mat_for_continuum_quantities_from_dofs_transpose(stuff.num_dofs, stuff.num_continuum_quantities);
try {calc_deriv_mats( 
    mat_for_continuum_quantities_from_dofs,
    mat_for_continuum_quantities_from_dofs_transpose,
    node_positions, 
    triangles, 
    stuff, 
    log);} catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}


// Calc masses assigned to dofs.
Eigen::Matrix<double,Eigen::Dynamic,1> dof_masses(3*stuff.num_nodes, 1);
try{calc_dof_masses(
    dof_masses,
    triangles, 
    stuff,
    log);} catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}


// Calculate some other things we'll need.
try {
    calc_char_mechanics_scales(triangles, stuff);
    calc_long_timescales(stuff, log);
    calc_damping_factor(stuff.damping_prefactor_1, stuff, log); // We'll recalculate later with stuff.damping_prefactor_2 when waiting for equil.
    calc_timestep(stuff, log);
    calc_steps_between_outputs(stuff, log);
    calc_dialling_procedure_params(stuff, log);
    } catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}


// This can be useful as a gentle way of allowing the simulation 
// to get going if all nodes in the initial data lay in a plane for
// example. Then, the random disturbance applied here will introduce small z
// components in the node positions, allowing the simulation to 'break out' of the
// plane. 
perturb_node_positions_with_noise( node_positions, stuff );

// Store the metric (a) and second fundamental form (b) of the reference state.
// In all these _comps matrices, each row contains the components xx, xy, yy.
Eigen::Matrix<double,Eigen::Dynamic,3> a_comps(stuff.num_tris, 3);
Eigen::Matrix<double,Eigen::Dynamic,3> b_comps(stuff.num_tris, 3);
a_comps.fill(1); 
a_comps.col(1).fill(0);
b_comps.fill(0);
Eigen::Matrix<double,Eigen::Dynamic,3> a_comps_at_zero_dial_factor = a_comps;
Eigen::Matrix<double,Eigen::Dynamic,3> b_comps_at_zero_dial_factor = b_comps;




// If an ansatz file was given, we move the node positions from those
// of the reference state to those given in the ansatz.
try{
    if( stuff.ansatz_file_was_given ){
        move_nodes_to_ansatz_positions(node_positions, ansatz_node_positions, stuff);
    }
} catch(const std::runtime_error &error){log << error.what() << std::endl; return -1;}


// Initialize dofs and many other things.
// matrix.block(i,j,p,q) means block of size (p,q) starting at (i,j).
dofs.block(0,0,stuff.num_nodes,1) = node_positions.col(0);
dofs.block(stuff.num_nodes,0,stuff.num_nodes,1) = node_positions.col(1);
dofs.block(2*stuff.num_nodes,0,stuff.num_nodes,1) = node_positions.col(2);
initialize_glass_slides(node_positions, stuff);
node_positions.resize(0,0);// Deleting node_positions matrix; henceforth the node positions are just stored in the dofs vector we just created.
Eigen::Matrix<double,Eigen::Dynamic,1> continuum_quantities(stuff.num_continuum_quantities, 1);
Eigen::Matrix<double,Eigen::Dynamic,3> normals(stuff.num_tris, 3);
Eigen::Matrix<double,Eigen::Dynamic,3> abar_comps(stuff.num_tris, 3);
Eigen::Matrix<double,Eigen::Dynamic,3> bbar_comps(stuff.num_tris, 3);
Eigen::Matrix<double,Eigen::Dynamic,1> forces(stuff.num_dofs, 1);
Eigen::Matrix<double,Eigen::Dynamic,1> velocities(stuff.num_dofs, 1);
Eigen::Matrix<double,Eigen::Dynamic,2> curvatures(stuff.num_tris, 2); // Gauss and mean.
Eigen::Matrix<double,Eigen::Dynamic,2> energy_densities(stuff.num_tris, 2); // Stretch and bend.
Eigen::Matrix<double,Eigen::Dynamic,1> angle_deficits(stuff.num_nodes, 1);
Eigen::Matrix<double,Eigen::Dynamic,8> stresses(stuff.num_tris, 8); // 8 elements in a row are eigenval_1, eigenvec_1 [3 comps], eigenval_2, eigenvec_2 [3 comps].
Eigen::Matrix<double,Eigen::Dynamic,1> strains(stuff.num_tris, 1);
std::vector<double> energies(3); // stretch_energy, bend_energy, kinetic_energy;
Eigen::Matrix<double,Eigen::Dynamic,1> def_shear_moduli(stuff.num_tris, 1);
Eigen::Matrix<double,Eigen::Dynamic,1> def_thicknesses(stuff.num_tris, 1);
Eigen::Matrix<double,Eigen::Dynamic,1> del_energy_by_del_continuum_quantities(stuff.num_continuum_quantities, 1);
velocities.fill(0.0); // Precautionary.
double time = 0.0;
int stepcount = 0;
double dial_factor = dial_factor_to_start_from;
int timesteps_since_previous_equil_check = -123456;
double last_dial_fac_equilibriated_at = dial_factor_to_start_from;
int timesteps_since_current_dial_phase_began = static_cast<int>( fmod(dial_factor_to_start_from,stuff.dial_factors_to_hold_at_spacing) * static_cast<double>(stuff.timesteps_per_dial_phase) ); // This is reset to zero every time equil is reached and a new `dialling' phase starts. Initializing it like this is cheeky, but if you stop and restart a simulation (using the last output of the first run as an ansatz) this ensures the amount of time until the next waiting_for_equil phase is as it was.
stuff.status = dialling; // Indicates whether the simulation is currently in a 'dialling' phase, or is not dialling and is instead waiting for equilibrium, or whether equilibrium has been reached but the next dialling phase has not yet begun.
// Set log format to something good for the regular outputs.
log << std::scientific << std::setprecision(8);
log << "Beginning dynamics loop." << std::endl;


// Begin dynamics loop.
while( true ){

    // Once the gradient descent dynamics (if any) have finished,
    // switch to Newtonian dynamics.
    if( stuff.using_gradient_descent_dynamics && stepcount > stuff.num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics ){
        log << "Finished doing gradient descent, now switching to Newtonian dynamics." << std::endl;
        stuff.using_gradient_descent_dynamics = false;
        // The following now need recalculating.
        calc_steps_between_outputs(stuff, log);
        calc_timestep(stuff, log);
    }
    

    // Dial_factor will end up very slightly
    // > 1 since it increases in discrete 
    // increments. Here we correct for that 
    // by replacing any >1 value with 1.
    if( dial_factor >= 1.0 ){
        dial_factor = 1.0;
    }


    // Switch from dialling to waiting-for-equil,
    // if it's time to switch.
    if( stuff.status == dialling && timesteps_since_current_dial_phase_began >= stuff.timesteps_per_dial_phase ){
        stuff.status = waiting_for_equil;
        // An EquilCheck has not actually just occurred, but this has
        // a desired effect of ensuring that waiting_for_equil phase 
        // is at least timesteps_between_equil_checks long. 
        // A hacky approach I know, and I feel bad about it.
        timesteps_since_previous_equil_check = 0.0;
        // Set damping factor to waiting-phase value.
        calc_damping_factor(stuff.damping_prefactor_2, stuff, log);
        calc_timestep(stuff, log); // Timestep depends on damping_factor for gradient_descent, so we recalculate (unlikely ever to be needed though).
        log << "Reached dial_factor = " << dial_factor << ", now waiting for equilibrium." << std::endl;
    }


    // Calculate geometry.
    try{
        calc_surface_derivs(
        continuum_quantities,
        mat_for_continuum_quantities_from_dofs,
        dofs);
    } catch(const std::runtime_error &error){log << "At stepcount = " << stepcount << ", dial_factor = " << dial_factor << " there was a problem: " << error.what() << std::endl; return -1;}
    
    

    // If dialling_from_ansatz_rather_than_ref_state, adjust 
    // the stored a_comps_at_zero_dial_factor and b_comps_at_zero_dial_factor 
    // accordingly.
    if( stepcount == 0 && stuff.dialling_from_ansatz_rather_than_ref_state ){
        calc_a_comps_and_b_comps_and_normals(
            a_comps, 
            b_comps,
            normals,
            continuum_quantities, 
            stuff);
        a_comps_at_zero_dial_factor = a_comps;
        b_comps_at_zero_dial_factor = b_comps;
    }

    // HERE WE DO THE DIALLING OF WHATEVER NEEDS DIALLING.
    // THIS IS PROBABLY THE PART OF THE CODE YOU'RE
    // MOST LIKELY TO WANT TO MODIFY.
    if( stuff.status == dialling || stepcount == 0){
        try{
            do_dialling(
            abar_comps,
            bbar_comps,
            def_thicknesses,
            def_shear_moduli,
            dial_factor, 
            abar_info, 
            bbar_info,
            a_comps_at_zero_dial_factor,
            b_comps_at_zero_dial_factor,
            triangles,
            nodes,
            dofs, 
            normals,
            log,
            stuff);}
        catch(const std::runtime_error &error){log << "At stepcount = " << stepcount << ", dial_factor = " << dial_factor << " there was a problem: " << error.what() << std::endl; return -1;}
    }


    // Calculate forces due to deformation.
    try{
    calc_deformation_forces(
        forces,
        del_energy_by_del_continuum_quantities,
        continuum_quantities,
        def_shear_moduli,
        def_thicknesses,
        abar_comps,
        bbar_comps,
        mat_for_continuum_quantities_from_dofs_transpose,
        triangles,
        stuff,
        log);}
        catch(const std::runtime_error &error){log << "At stepcount = " << stepcount << ", dial_factor = " << dial_factor << " there was a problem: " << error.what() << std::endl; return -1;}
    


    // Calculate non-deformation forces, e.g. 
    // damping, weight.
    try{
    calc_non_deformation_forces(
        forces,
        velocities,
        dofs,
        dof_masses,
        stuff);}
        catch(const std::runtime_error &error){log << "At stepcount = " << stepcount << ", dial_factor = " << dial_factor << " there was a problem: " << error.what() << std::endl; return -1;}


    // Impose any constraints, e.g. clamped boundaries.
    impose_constraints(
        nodes,
        dofs,
        velocities,
        forces,
        triangles,
        stuff);

    // Check for equilibrium
    if( stuff.status == waiting_for_equil && timesteps_since_previous_equil_check > stuff.timesteps_between_equil_checks ){
        log << "Checking for equilibrium at " << get_real_time() << ", stepcount = " << stepcount << ", simulation time = " << time << ", dial_factor = " << dial_factor << std::endl;
        if( is_at_equilibrium(
                forces,
                velocities,
                stuff,
                log) ){
            stuff.status = equil_reached;
            log << "Equilibrium reached. Writing VTK output at " << get_real_time() << ", stepcount = " << stepcount << ", simulation time = " << time << ", dial_factor = " << dial_factor << std::endl;    
            last_dial_fac_equilibriated_at = dial_factor; }
        timesteps_since_previous_equil_check = 0;
    }
    
    // Write data output at regular
    // intervals, and also whenever
    // equil is reached.
    if( (stepcount % stuff.steps_between_outputs == 0 && stuff.output_regularly) || stuff.status == equil_reached ){
            calc_a_comps_and_b_comps_and_normals(
                a_comps, 
                b_comps,
                normals,
                continuum_quantities, 
                stuff);
        calc_curvatures(curvatures, continuum_quantities, b_comps, stuff);
        calc_angle_deficits(angle_deficits, dofs, nodes, triangles, stuff);
        calc_energy_densities(energy_densities, a_comps, b_comps, abar_comps, bbar_comps, def_shear_moduli, def_thicknesses, stuff, log);
        calc_energies(energies, energy_densities, velocities, dof_masses, triangles, stuff);
        calc_stresses(stresses, continuum_quantities, abar_comps, bbar_comps, def_shear_moduli, def_thicknesses, triangles, stuff, log);
        calc_strains(strains, continuum_quantities, a_comps, abar_comps, stuff);
        try{write_output(
            dofs, 
            triangles, 
            stepcount, 
            dial_factor, 
            curvatures,
            energy_densities, 
            angle_deficits, 
            stresses,
            strains,
            energies,       
            stuff);} catch(const std::runtime_error &error){log << "At stepcount = " << stepcount << ", dial_factor = " << dial_factor << " there was a problem: " << error.what() << std::endl; return -1;}
        log << "Wrote VTK output at " << get_real_time() << ", stepcount = " << stepcount << ", simulation time = " << time+stuff.timestep << ", current dial factor = " << dial_factor << std::endl;
    }


    // Advance dynamics (i.e. time integration).
    try{advance_dynamics(
            dofs, 
            velocities,
            forces,
            dof_masses,
            stuff,
            log);} catch(const std::runtime_error &error){
                calc_a_comps_and_b_comps_and_normals(
                    a_comps, 
                    b_comps,
                    normals,
                    continuum_quantities, 
                    stuff);
                calc_curvatures(curvatures, continuum_quantities, b_comps, stuff);
                calc_angle_deficits(angle_deficits, dofs, nodes, triangles, stuff);
                calc_energy_densities(energy_densities, a_comps, b_comps, abar_comps, bbar_comps, def_shear_moduli, def_thicknesses, stuff, log);
                calc_energies(energies, energy_densities, velocities, dof_masses, triangles, stuff);
                calc_stresses(stresses, continuum_quantities, abar_comps, bbar_comps, def_shear_moduli, def_thicknesses, triangles, stuff, log);
                calc_strains(strains, continuum_quantities, a_comps, abar_comps, stuff);
                write_output(
                    dofs, 
                    triangles, 
                    stepcount, 
                    dial_factor, 
                    curvatures,
                    energy_densities, 
                    angle_deficits, 
                    stresses,
                    strains,
                    energies,       
                    stuff);
                    log << "At stepcount = " << stepcount << ", dial_factor = " << dial_factor << " there was a problem: " << error.what() << std::endl; return -1;
                }



    time += stuff.timestep;
    stepcount += 1;
    if( stuff.status == dialling){
        timesteps_since_current_dial_phase_began += 1;
    }
    dial_factor = last_dial_fac_equilibriated_at 
                    + stuff.dial_factors_to_hold_at_spacing * 
                        (static_cast<double>(timesteps_since_current_dial_phase_began) / static_cast<double>(stuff.timesteps_per_dial_phase));
    if( stuff.status == waiting_for_equil){
        timesteps_since_previous_equil_check += 1;
    }


    // If equilibrium has been reached and we are at the end of the 
    // simulation, terminate the while(true) loop.
    if( dial_factor >= 1.0 && stuff.status == equil_reached && stuff.slide_stiffness_prefactor > 0.0 ){ 
        break;
    }


    // If equil has been reached but dial_factor hasn't
    // reached 1 yet, move on to the next dialling
    // phase.
    if( stuff.status == equil_reached ){
        stuff.status = dialling;
        timesteps_since_current_dial_phase_began = 0;
        // Set damping factor to dialling-phase value.
        calc_damping_factor(stuff.damping_prefactor_1, stuff, log);
        calc_timestep(stuff, log); // Timestep depends on damping_factor for gradient_descent, so we recalculate (unlikely ever to be needed though).
    }
}



// Print some helpful final things.
std::cout << "\n \n " << std::endl;
log << "Simulation completed at time = " << time << ", stepcount = " << stepcount << std::endl;

return 0;    
}