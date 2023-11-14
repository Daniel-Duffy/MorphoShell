/* 
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>
#include <vector>

#include "Node_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"



bool is_at_equilibrium(
    const Eigen::Matrix<double,Eigen::Dynamic,1> &forces,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Stuff_Class &stuff,
    Out_Stream_Class &log);