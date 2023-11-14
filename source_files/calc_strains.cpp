/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/



//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>
#include <omp.h>


#include "calc_strains.hpp"
#include "Stuff_Class.hpp"

void calc_strains(
    Eigen::Matrix<double,Eigen::Dynamic,1> &strains,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &a_comps,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &abar_comps,
    const Stuff_Class &stuff){

    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // Index into continuum_quantities that is the start of triangle t's set of 15 continuum quantities
        int q = 15*t;

        double Xx  = continuum_quantities(q+0);
        double Xy  = continuum_quantities(q+1);
        double Yx  = continuum_quantities(q+5);
        double Yy  = continuum_quantities(q+6);
        double Zx  = continuum_quantities(q+10);
        double Zy  = continuum_quantities(q+11);

        Eigen::Matrix<double,3,2> def_gradient{{Xx, Xy},
                                               {Yx, Yy},
                                               {Zx, Zy}};

        // This could also be calculated as def_gradient.transpose()*def_gradient.
        Eigen::Matrix<double,2,2> a {{a_comps(t,0), a_comps(t,1)},
                                    {a_comps(t,1), a_comps(t,2)},};

        Eigen::Matrix<double,2,2> abar {{abar_comps(t,0), abar_comps(t,1)},
                                        {abar_comps(t,1), abar_comps(t,2)},};


        Eigen::Matrix<double,2,3> def_gradient_pseudo_inv = (def_gradient.transpose() * def_gradient).inverse() * (def_gradient.transpose());

        // This is a tensor that you can dot on either side with a unit-length
        // direction IN THE TANGENT PLANE, and it will spit our a scalar that
        // to leading order equals the strain relative to the programmed metric,
        // i.e. (dx - d\bar{X}) / d\bar{X} where dx is a deformed-state
        // infinitesimal length in the direction you chose, and d\bar{X} is the
        // `desired'/`programmed` length of this vector according to the
        // progMetric. See my notes on A Strain Tensor For Active Shells.
        Eigen::Matrix<double,3,3> my_strain_tensor = 0.5 * (def_gradient_pseudo_inv.transpose() * (a - abar) * def_gradient_pseudo_inv);

        // To boil the characteristic size of this tensor down into a scalar we
        // need to apply some kind of matrix norm. We choose to use 1/sqrt(2) * the
        // Frobenius norm, because that corresponds to the RMS of the two
        // principal strains (eigenvalues of our strain tensor).
        strains(t) = sqrt( 0.5 * ((my_strain_tensor.transpose() * my_strain_tensor).trace()) );


    }
}
