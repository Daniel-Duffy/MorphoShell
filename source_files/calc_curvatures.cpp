/////////////////////////////////////////////////////
/*
Copyright (C) 2023, Daniel Duffy, daniellouisduffy@gmail.com. All rights reserved.
Please cite Daniel Duffy and John S. Biggins if you 
use any part of this code in work that you publish or distribute.

This file is part of MorphoShell.

MorphoShell is distributed under the terms of the Cambridge Academic
Software License (CASL). You should have received a copy of the license
along with MorphoShell. If not, contact Daniel Duffy, daniellouisduffy@gmail.com.
*/
/////////////////////////////////////////////////////

// Turn Eigen bounds checking off for speed.
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "calc_curvatures.hpp"
#include "Node_Class.hpp"
#include "Tri_Class.hpp"
#include "Stuff_Class.hpp"
#include "kahan_sum.hpp"

void calc_curvatures(
    Eigen::Matrix<double,Eigen::Dynamic,2> &curvatures,
    const Eigen::Matrix<double,Eigen::Dynamic,1> &continuum_quantities,
    const Eigen::Matrix<double,Eigen::Dynamic,3> &b_comps,
    const Stuff_Class &stuff){

    #pragma omp parallel for simd
    for(int t = 0; t < stuff.num_tris; ++t){

        // Index into continuum_quantities that is the start of triangle t's set of 15 continuum quantities
        int q = 15*t;

        // continuum_quantities structure will be 
        // (tri1_Xx, tri1_Xy, tri1_Xxx, tri1_Xxy, tri1_Xyy, tri1_Yx, tri1_Yy, tri1_Yxx, tri1_Yxy, tri1_Yyy, tri1_Zx, tri1_Zy, tri1_Zxx, tri1_Zxy, tri1_Zyy, tri2_Xx, ...)
        //     0        1         2        3          4       5         6        7         8         9        10       11        12       13         14
        double Xx  = continuum_quantities(q+0);
        double Xy  = continuum_quantities(q+1);
        double Yx  = continuum_quantities(q+5);
        double Yy  = continuum_quantities(q+6);
        double Zx  = continuum_quantities(q+10);
        double Zy  = continuum_quantities(q+11);

        // A direct approach would just be to calc the shape operator, which
        // is a.inverse() * b, and then take its det and trace; the former 
        // being the Gauss curvature, and the latter being twice the mean
        // curvature.
        // Instead we use the 2x2 'symmetric shape operator', which should have
        // better numerical properties - see Wang, Clark, Jiao 2009 and 
        // Jiao, Zha 2008. It works by essentially factoring out and thereby 
        // ignoring the pure rotation part of the deformation, which has no 
        // effect on curvatures. This is accomplished via QR factorisation 
        // of the deformation gradient following the 2009 paper, though the 
        // 2008 one uses SVD instead.

        Eigen::Matrix<double,3,2> def_gradient{{Xx, Xy},
                                               {Yx, Yy},
                                               {Zx, Zy},};

        Eigen::Matrix<double,2,2> secFF {{b_comps(t,0), b_comps(t,1)},
                                         {b_comps(t,1), b_comps(t,2)},};        

        Eigen::HouseholderQR< Eigen::Matrix<double,3,2> > qr_decomp(3, 2);
        qr_decomp.compute(def_gradient);

        // Extracting the 'R' part relies on the special and under-documented
        // way Eigen stores the QR-decomposed matrix, which I believe is the LAPACK
        // way. Inconveniently, it returns R as 3x2 upper triangular, not 2x2.
        Eigen::Matrix<double,3,2> stretchPart_3x2 = qr_decomp.matrixQR().triangularView<Eigen::Upper>();
        Eigen::Matrix<double,2,2> invStretchPart = stretchPart_3x2.block<2,2>(0,0).inverse(); // If it's not invertible you'll notice elsewhere!

        Eigen::Matrix<double,2,2> sym_shape_op = invStretchPart.transpose() * secFF * invStretchPart;

        curvatures(t,0) = sym_shape_op.determinant(); // Gauss
        curvatures(t,1) = 0.5 * sym_shape_op.trace(); // Mean
        // In case you want principle curvatures too, you could follow eqn 7 in the
        // 2009 paper --- see the commented-out lines just below.
        /*
        double sqrt_discrim = sqrt( (sym_shape_op(0,0) - sym_shape_op(1,1)) * (sym_shape_op(0,0) - sym_shape_op(1,1)) + 4.0 * sym_shape_op(0,1) * sym_shape_op(0,1) );
        double princ_curv1 = 0.5 * (sym_shape_op(0,0) + sym_shape_op(1,1) + sqrt_discrim);
        double princ_curv2 = 0.5 * (sym_shape_op(0,0) + sym_shape_op(1,1) - sqrt_discrim);
        */
        
        // For principal curvature directions, you could follow eqn(8) in the 
        // 2009 paper (in which the \hat{u} vecs should be \hat{q} vecs --- see the
        // bottom of their table 1).
    }
}
