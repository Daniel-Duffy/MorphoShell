
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "Out_Stream_Class.hpp"


void calc_deriv_mats(
    Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_for_continuum_quantities_from_dofs,
    Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_for_continuum_quantities_from_dofs_transpose,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions,
    const std::vector<Tri_Class> &triangles,
    const Stuff_Class &stuff,
    Out_Stream_Class &log);