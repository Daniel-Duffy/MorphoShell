
# Makefile by Daniel Duffy, <daniellouisduffy@gmail.com>
# For help, email me or see: 
# https://www.gnu.org/software/make/manual/html_node/Prerequisite-Types.html#Prerequisite-Types

# To include headers not in the same directory as your .cpp files, uncomment and adjust:
# INCDIRS = -I./someDirOfHeaders -I/usr/local/someOtherDirOfHeaders
# To view the default include paths, I found the following terminal 
# command on Stackoverflow 11946294
# g++-14  -E -x c++ - -v < /dev/null

# C Pre-processor flags
CPPFLAGS = $(INCDIRS)

# C++ compiler flags
# The -march=native flag may need to removed/replaced to compile on Windows.
CXXFLAGS = -fopenmp -std=c++17 -pedantic -g -march=native -O3\
-Wall -Wconversion -Wextra\
-Wformat  -Wshadow -Werror\
-Wpointer-arith -Wcast-qual -Wwrite-strings\

# If you need to link to particular libraries, not in the default search path. This can be 
# required if you built and installed gcc manually for example, which 
# for me was installed to /opt/gcc-9.3.0
# -Wl,-rpath=/path/to/mylibdir adds mylibdir to the start of the LINK-TIME library search path.
# -Lmylibdir adds mylibdir to the start of the RUN-TIME library search path.
# LIBS = -lconfig++ -Wl,-rpath=/opt/gcc-9.3.0/lib64 -L/opt/gcc-9.3.0/lib64
# If you don't need the above, e.g. you installed gcc normally via sudo apt install, this should be fine:
# LIBS = -lconfig++ 
LIBS = -fopenmp

# This specifies the compiler command (executable) to execute. So for a normal sudo apt install gcc-14
# installation this should be fine:
CXX = g++-14
# But if you had to install gcc-9 manually you *may* need something like this instead:
# CXX = /opt/gcc-9.3.0/bin/g++

# Path to source files subdirectory relative to directory holding Makefile.
VPATH = ./source_files

# Name of executable file to make.
EXEFILE = morphoshell.exe


# Directory to put .o files in to avoid mess.
OBJDIR := object_files

# List of object files which will be created (functions, classes etc).
# Order should not matter on most modern compilers, but good practice 
# is to have fnA.o before fnB.o if fnA calls fnB, as with libraries where
# order matters more often. For discussion see:
# An Introduction to GCC, Brian Gough.
OBJECTS := $(addprefix $(OBJDIR)/,\
main.o \
get_real_time.o \
extract_filename.o \
read_settings.o \
Out_Stream_Class.o \
Stuff_Class.o \
initialize_directories.o \
preliminary_setup.o \
read_vtk_data.o \
calc_mesh_info.o \
calc_triangles_incident_on_nodes.o \
calc_node_neighbours.o \
calc_triangle_adjacencies.o \
calc_edges.o \
Node_Class.o \
Tri_Class.o \
Edge_Class.o \
kahan_sum.o \
set_up_patch_fitting.o \
calc_deriv_mats.o \
calc_long_timescales.o \
calc_damping_factor.o \
calc_timestep.o \
calc_steps_between_outputs.o \
calc_char_mechanics_scales.o \
calc_dialling_procedure_params.o \
perturb_node_positions_with_noise.o \
calc_dof_masses.o \
move_nodes_to_ansatz_positions.o \
initialize_glass_slides.o\
calc_surface_derivs.o \
calc_a_comps_and_b_comps_and_normals.o \
do_dialling.o \
calc_deformation_forces.o \
calc_non_deformation_forces.o \
find_extremum_of_doubles.o \
is_at_equilibrium.o \
impose_constraints.o \
calc_curvatures.o \
calc_angle_deficits.o \
calc_energy_densities.o \
calc_energies.o \
calc_stresses.o \
calc_strains.o \
write_output.o \
advance_dynamics.o \
)

# For big compile jobs, you can list all header file dependencies for each object file in turn, 
# using commands like g++ -MM main.o to handle dependencies to avoid recompiling the whole thing every
# time. This sounds complicated to me, but is certainly doable, see here:
# https://www.gnu.org/software/make/manual/html_node/Automatic-Prerequisites.html#Automatic-Prerequisites
# Certainly, not doing this can lead to problems pretty quickly (and has for me).
# I recommend what I see as the easy solution for a small code like this:
# just run "make clean" every time before you recompile, to start from scratch 
# each time, unless you've only made a small change in a single .cpp file (not a header).

.PHONY: clean
clean:
#	rm *.o  $(EXEFILE) # from older version
	rm -r $(OBJDIR) $(EXEFILE)
	
	
$(OBJECTS): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(EXEFILE) : $(OBJECTS)
	${CXX} $^ $(LIBS) -o $@
