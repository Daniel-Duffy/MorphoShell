//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "kahan_sum.hpp"

void calc_energies(
    std::vector<double> &energies,
    const Eigen::Matrix<double,Eigen::Dynamic,2> &energy_densities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &velocities,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dof_masses,
    const std::vector< Tri_Class > &triangles,
    const Stuff_Class &stuff);