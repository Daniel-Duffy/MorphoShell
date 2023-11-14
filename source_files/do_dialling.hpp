
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <cstddef>

#include "Tri_Class.hpp"
#include "Node_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"

void do_dialling(
    Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    Eigen::Matrix<double,Eigen::Dynamic,3> &bbar_comps,
    Eigen::Matrix<double,Eigen::Dynamic,1> &def_thicknesses,
    Eigen::Matrix<double,Eigen::Dynamic,1> &def_shear_moduli,
    double dial_factor,
    const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &abar_info,
    [[maybe_unused]] const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &bbar_info,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &initial_a_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &initial_b_comps,
    [[maybe_unused]] const std::vector<Tri_Class> &triangles,
    [[maybe_unused]] const std::vector<Node_Class> &nodes,
    [[maybe_unused]] const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs,
    [[maybe_unused]] const Eigen::Matrix<double,Eigen::Dynamic,3> &normals,
    [[maybe_unused]] Out_Stream_Class &log,
    const Stuff_Class &stuff);