/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/


#include <iostream>
#include <cmath>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void advance_dynamics(
    Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    Stuff_Class &stuff,
    Out_Stream_Class &log);