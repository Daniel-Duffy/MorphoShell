
//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>

void calc_surface_derivs(
    Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat_for_continuum_quantities_from_dofs,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &dofs);