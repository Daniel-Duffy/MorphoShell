
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "Stuff_Class.hpp"

void calc_energy_densities(
    Eigen::Matrix<double,Eigen::Dynamic,2> &energy_densities,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &a_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &b_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &bbar_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &def_shear_moduli,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &def_thicknesses,
    const Stuff_Class &stuff);