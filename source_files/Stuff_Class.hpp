/*
/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

This is the header file for the 'stuff' class that will hold many useful things. 
It should always be passed to functions by reference, and usually const.*/

#ifndef _STUFF_CLASS_TAG_
#define _STUFF_CLASS_TAG_

#include <string> // Used for creating directory for output data

#include "Status_Enum.hpp"

class Stuff_Class{

public:

	// SETTINGS, WHICH WILL BE READ
	// IN FROM THE SETTINGS FILE.
	double outputs_per_char_long_time;
	bool dialling_from_ansatz_rather_than_ref_state;
	double patch_mat_dimless_conditioning_thresh;
	double ref_thickness_if_uniform;
	double def_poisson_ratio; // Deformed-state!
	double ref_shear_modulus_if_uniform; // Deformed-state!
	double ref_density;
	double dial_phase_time_prefactor;
	double time_between_equil_checks_prefactor;
	double damping_prefactor_1;
	double damping_prefactor_2;
	double timestep_prefactor;
	double dial_factors_to_hold_at_spacing;
	double force_to_char_force_ratio_equil_threshold;
	double speed_to_char_speed_ratio_equil_threshold;
	double bend_stiffening_scale_factor;
	double slide_stiffness_prefactor;
	double slide_friction_coefficent;
	double slide_speed_prefactor;
	bool using_initial_glass_slide_z_coords_from_settings_file;
	double initial_lower_slide_z_coord;
	double initial_upper_slide_z_coord;


	// OTHER USEFUL QUANTITIES THAT
	// ARE NOT READ IN, BUT RATHER 
	// ARE CALCULATED BY MORPHOSHELL
	// ITSELF AT RUN-TIME.
	std::string settings_filepath;
	std::string ref_data_filepath;
	std::string ansatz_filepath;
	std::string output_dir_path;
	std::string aux_output_filepath;
	bool ansatz_file_was_given;

	bool output_regularly;

	double sample_char_length;
	double approx_min_tri_size;

	double min_ref_thickness;
	double min_ref_shear_modulus;
	double stretch_long_time;
	double bend_long_time;
	double char_long_time;
	double dial_phase_time;
	long int timesteps_between_equil_checks;
	double damping_factor;
	double timestep;
	long int steps_between_outputs;
	double char_flexural_rigidity;
	double char_force_scale;
	double char_speed_scale;
	double char_stretch_energy_density_scale;
	double char_stretch_energy_scale;
	double total_ref_area;

	double upper_slide_z_coord;
	double upper_slide_tot_vert_force;
	double lower_slide_z_coord;
	double lower_slide_tot_vert_force;

	double grav_field_strength;

	int num_steps_to_use_gradient_descent_dynamics_for_before_switching_to_newtonian_dynamics; 
	bool using_gradient_descent_dynamics;


	//double dial_factor_change_per_timestep;
	long int timesteps_per_dial_phase;

	Status_Enum status;

    // Total numbers of nodes, triangles, and edges, etc
    int num_nodes, num_boundary_nodes, num_tris, num_boundary_tris, num_edges, num_boundary_edges, num_dofs, num_continuum_quantities;
};
#endif
