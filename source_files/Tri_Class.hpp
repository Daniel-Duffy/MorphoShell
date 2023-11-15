/*
/////////////////////////////////////////////////////

/////////////////////////////////////////////////////

This is the header file for the class that will contain data for each
triangular element, such as vertices, area etc.*/

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

#ifndef _TRI_CLASS_TAG_
#define _TRI_CLASS_TAG_

#include <vector>
#include <Eigen/Dense>

#include "Out_Stream_Class.hpp"
#include "Node_Class.hpp"

class Tri_Class{
public:

    double ref_thickness;
    double ref_shear_modulus;

    //double normal_sign_fac;

    // Current unit normal to face.
    Eigen::Vector3d normal;

    double ref_area;
    double area;

    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
    Out_Stream_Class tri_log_stream;

    // id so this triangle 'knows' which it is.
    int id;

    // tag that will be read in from input file. Makes it easy to mark out certain
    // triangles for special treatment.
    bool tag;

    /* Boolean representing whether the triangle has a vertex on the boundary of
    the sample (true) or not (false).*/
    bool is_boundary;

    // Initial (reference) area and 1/(current area) for this triangle.
    double initArea;
    double invCurrArea;

    /* Indices (ids) of the nodes at the vertices of the triangle. These are
    in no particular order*/
    std::array<int, 3> vertex_ids;

    /* ids (and indices in the triangles container vector) of the 3 edges of
    this triangle.*/
    std::array<int, 3> edge_ids;

    /* Vector storing, for each non-boundary edge, its initial
    length divided by sum of initial non-boundary edge lengths for the triangle.
    These fractions are used for weighting the edge normals used in the secFF
    calculation.*/
    Eigen::VectorXd initNonBoundEdgeLengthFracs;

    /* ids (and indices in the triangles container vector) of the other
    triangles that share an edge with this triangle.*/
    std::vector<int> edge_sharing_tri_ids;

    /* A vector with one component for each non-boundary edge of this triangle,
    with each entry being either +1.0 or -1.0. A value of +1.0 means that
    this triangle corresponds to adjTriids(0) for the corresponding edge.
    -1.0 implies adjTriids(1) similarly. */
    Eigen::VectorXi edgeAdjTriidSelectors;

    /* Vector where the i'th element is the index corresponding to this triangle
    in the edge_sharing_tri_ids vector of the i'th edge-sharing triangle
    (neighbour) of THIS triangle.*/
    Eigen::VectorXi indicesIntoEdgeSharingTriidsOfNeighbours;

    /* Indices (ids) of the 3 nodes that are not vertices of the triangle, but
    are part of the estimation of the 2nd F.F. for each triangle. For
    non-boundary triangles these are the non-shared-edge nodes from each of the
    triangles sharing an edge with this triangle. For boundary triangles one or
    two of these are missing, and different nodes are chosen based on proximity
    to this triangle's centroid. There is again no particular order.*/
    std::array<int, 3>  non_vertex_patch_nodes_ids;

    /* Position of triangle's centroid in initial (reference) x-y plane.*/
    Eigen::Vector3d ref_centroid;

    /* Matrix (of doubles) where each *column* is a vector describing a side of
    the triangle. We only need two sides per triangle for the algorithm
    hence the 3x2 matrices.*/
    Eigen::Matrix<double,3,2> currSides;

    // Current unit normal to face.
    Eigen::Vector3d faceNormal;

    /* Inverse of 2x2 matrix that has (two) initial sides of the triangle as
    columns. Those sides correspond to the two current sides stored in currSides.*/
    Eigen::Matrix<double,2,2> invInitSidesMat;

    /* The reference (initial) state in-triangle-plane *outward* normals of
    each triangle's sides (with lengths = corresponding side lengths).
    initOutwardTriNormals.col(v) is opposite vertex_ids(v).*/
    Eigen::Matrix<double,2,3> initOutwardSideNormals;

    /* Matrix representing the inverse of the *energetically* favoured metric for
    this triangle, induced by a programmed nematic director field, for example.
    The matrix stored here is the 'dialled in' value, which is used to prevent
    anything too explosive happening. */
    Eigen::Matrix<double,2,2> dialledInvProgMetric;

    // Determinant of the above dialledInvProgMetric matrix.
    double detDialledInvProgMetric;

    /* Dialled in programmed scalar 'tau' factor.*/
    double dialledProgTau;

    /* Matrix representing the components of the *energetically* favoured
    (programmed) Second Fundamental Form in the x-y cartesian coordinate system
    of the initial flat state. */
    Eigen::Matrix<double,2,2> dialledProgSecFF;

    /* Matrix representing the total deformation gradient of this triangle. This
    maps the in-plane reference (initial) state (for which no z components are
    stored) to a triangle in 3D space, so it is 3x2.*/
    Eigen::Matrix<double,3,2> defGradient;

    /* The (1st order) approximation for the metric for this triangle, which
    equals defGradient.transpose9) * defGradient.*/
    Eigen::Matrix<double,2,2> metric;

    /* Inverse of metric.*/
    Eigen::Matrix<double,2,2> invMetric;

    /* Determinant of the inverse of the metric.*/
    double detInvMetric;

    /* Matrix representing (1st Piola-Kirchoff stress tensor)/2 for this triangle.*/
    Eigen::Matrix<double,3,2> halfPK1Stress;

    /*Matrix that is pre-calculated and then used repeatedly in finding the
    components of the second fundamental form estimated for this triangle. */
    Eigen::Matrix<double,6,3> mat_for_patch_2nd_derivs;

    /* Matrix of the second position derivatives, used in calculating the
    secFF estimate. */
    Eigen::Matrix<double,3,3> patchSecDerivs;

    /* Estimated Second Fundamental Form matrix (secFF) of the deformed surface,
    defined (as with the deformation gradient) with respect to the 'material'
    coodinate system, i.e. the coordinate chart that used to be the (x,y)
    cartesians of the flat initial state, and then deformed with the sheet. This
    is a 2x2 symmetric matrix. */
    Eigen::Matrix<double,2,2> secFF;

    /* Derivative of the bending energy density with respect to the secFF.*/
    Eigen::Matrix<double,2,2> energyDensityDerivWRTSecFF;

    /* The derivative of the bending energy density with respect to the metric.*/
    Eigen::Matrix<double,2,2> bendEnergyDensityDerivWRTMetric;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for log).
    I ensure that default data values are recognisable values,
    for debugging.
    edge_sharing_tri_ids, nonSharedVerticesOfEdgeSharingTris,
    edgeAdjTriidSelectors, indicesIntoEdgeSharingTriidsOfNeighbours,
    and usefulTermsForSecFFDeriv
    are left with zero size at initialisation. */
    /*Tri_Class(){
        id = -12345;
        is_boundary = false;
        //vertex_ids.fill(123456789);
        non_vertex_patch_nodes_ids.fill(987654321);
        //edge_ids.fill(-4321);
        initArea = 98765;
        invCurrArea = 56789;
        currSides.fill(5678);
        faceNormal.fill(8765);
        initOutwardSideNormals.fill(8765432);
        invInitSidesMat.fill(-5678);
        dialledInvProgMetric.fill(-6789);
        detDialledInvProgMetric = -9.87654321;
        dialledProgTau = 123456;
        dialledProgSecFF.fill(-5678);
        defGradient.fill(123456);
        metric.fill(1234567.89);
        detInvMetric = -1234.567;
        invMetric.fill(-1234567.89);
        halfPK1Stress.fill(-1234567);
        patchSecDerivs.fill(12345.6789);
        secFF.fill(-7654);
        energyDensityDerivWRTSecFF.fill(456789);
        bendEnergyDensityDerivWRTMetric.fill(-456789);
        matForPatchSecDerivs.fill(54321);
    }
    */

    // Declare other member functions.

    // Calculate area.
    void calc_area(const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions);

    // Debugging function to display all member data.
    void display();

};
#endif
