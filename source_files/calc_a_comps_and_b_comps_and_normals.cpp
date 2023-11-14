/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

*/

//Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>
#include <omp.h>


#include "calc_a_comps_and_b_comps_and_normals.hpp"
#include "Stuff_Class.hpp"

void calc_a_comps_and_b_comps_and_normals(
    Eigen::Matrix<double,Eigen::Dynamic,3> &a_comps,
    Eigen::Matrix<double,Eigen::Dynamic,3> &b_comps,
    Eigen::Matrix<double,Eigen::Dynamic,3> &normals,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Stuff_Class &stuff){

    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // Index into continuum_quantities that is the start of triangle t's set of 15 continuum quantities
        int q = 15*t;

        // continuum_quantities structure will be (tri1_Xx, tri1_Xy, tri1_Xxx, tri1_Xxy, tri1_Xyy, tri1_Yx, tri1_Yy, tri1_Yxx, tri1_Yxy, tri1_Yyy, tri1_Zx, tri1_Zy, tri1_Zxx, tri1_Zxy, tri1_Zyy, tri2_Xx, ...)
        Eigen::Vector3d R_x  = {continuum_quantities(q+0), continuum_quantities(q+5), continuum_quantities(q+10)};
        Eigen::Vector3d R_y  = {continuum_quantities(q+1), continuum_quantities(q+6), continuum_quantities(q+11)};
        Eigen::Vector3d R_xx = {continuum_quantities(q+2), continuum_quantities(q+7), continuum_quantities(q+12)};
        Eigen::Vector3d R_xy = {continuum_quantities(q+3), continuum_quantities(q+8), continuum_quantities(q+13)};
        Eigen::Vector3d R_yy = {continuum_quantities(q+4), continuum_quantities(q+9), continuum_quantities(q+14)};


        // Metric with both indices down.
        a_comps(t,0) = R_x.dot(R_x);
        a_comps(t,1) = R_x.dot(R_y);
        a_comps(t,2) = R_y.dot(R_y);

        // Surface normal vector.
        Eigen::Vector3d temp = R_x.cross(R_y);
        normals.row(t) = temp / temp.norm();

        // Curvature tensor with both indices down, i.e. second fundamental form.
        b_comps(t,0) = -normals.row(t).dot(R_xx);
        b_comps(t,1) = -normals.row(t).dot(R_xy);
        b_comps(t,2) = -normals.row(t).dot(R_yy);
    }
}
