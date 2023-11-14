/* 
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this 
code in work that you publish or distribute.

This file is part of Shellmorph.

Shellmorph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Shellmorph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

This file defines the member functions for the Tri_Class class that holds
the data for each triangular element. The class constructors are the exception;
they are left in the header file for clarity there.*/

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include <vector>

#include "Tri_Class.hpp"
#include "Node_Class.hpp"

// Calculate area.
void Tri_Class::calc_area(const Eigen::Matrix<double,Eigen::Dynamic,3> &node_positions){
    Eigen::Vector3d side1 = node_positions.row(vertex_ids[1]) - node_positions.row(vertex_ids[0]);
    Eigen::Vector3d side2 = node_positions.row(vertex_ids[2]) - node_positions.row(vertex_ids[0]);
    area = (0.5 * side1.cross(side2)).norm();
}


//This is a debugging tool to display the node's data
/*
void Tri_Class::display() {
    tri_log_stream.open();
    tri_log_stream << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    tri_log_stream << "Triangle " << id << ":" << std::endl;
    tri_log_stream << "Boundary indicator = " << is_boundary << std::endl;
    tri_log_stream << "Initial (reference) area = " << initArea << std::endl;
    tri_log_stream << "invCurrArea = " << invCurrArea << std::endl;
    //tri_log_stream << "ids of vertices: " << vertex_ids.transpose() << std::endl;
    //tri_log_stream << "ids of edges: " << edge_ids.transpose() << std::endl;
    tri_log_stream << "Initial non-boundary edge length fractions: " << initNonBoundEdgeLengthFracs.transpose() << std::endl;
    //tri_log_stream << "ids of adjacent (edge-sharing) triangles: " << edge_sharing_tri_ids.transpose() << std::endl;
    //tri_log_stream << "edgeAdjTriidSelectors: " << edgeAdjTriidSelectors.transpose() << std::endl;
    tri_log_stream << "indicesIntoEdgeSharingTriidsOfNeighbours: " << indicesIntoEdgeSharingTriidsOfNeighbours.transpose() << std::endl;
    //tri_log_stream << "ids of non-vertex nodes in this triangle's patch: " << nonVertexPatchNodesids.transpose() << std::endl;
    tri_log_stream << "Current sides = " << "\n" << currSides << std::endl;
    tri_log_stream << "faceNormal = " << "\n" << faceNormal << std::endl;
    tri_log_stream << "Initial outward side normals = " << "\n" << initOutwardSideNormals << std::endl;
    tri_log_stream << "invInitInPlaneSidesMat = " << "\n" << invInitSidesMat << std::endl;
    tri_log_stream << "Dialled in inverse of programmed metric = " << "\n" << dialledInvProgMetric << std::endl;
    tri_log_stream << "detDialledInvProgMetric = " << "\n" << detDialledInvProgMetric << std::endl;
    tri_log_stream << "Dialled in prog tau factor = " << dialledProgTau << std::endl;
    tri_log_stream << "Dialled in programmed second fundamental form = " << "\n" << dialledProgSecFF << std::endl;
    tri_log_stream << "Deformation gradient = " << "\n" << defGradient << std::endl;
    tri_log_stream << "metric = " << "\n" << metric << std::endl;
    tri_log_stream << "Inverse of Metric = " << "\n" << invMetric << std::endl;
    tri_log_stream << "Det of inverse of Metric = " << "\n" << detInvMetric << std::endl;
    tri_log_stream << "matForPatchSecDerivs = " << "\n" << matForPatchSecDerivs << std::endl;
    tri_log_stream << "patchSecDerivs = " << "\n" << patchSecDerivs << std::endl;
    tri_log_stream << "Second fundamental form = " << "\n" << secFF << std::endl;
    tri_log_stream << "Bending energy density deriv wrt secFF = " << "\n" << energyDensityDerivWRTSecFF << std::endl;
    tri_log_stream << "bendEnergyDensityDerivWRTMetric = " << "\n" << bendEnergyDensityDerivWRTMetric << std::endl;
    tri_log_stream << "halfPK1Stress = " << "\n" << halfPK1Stress << std::endl;
    tri_log_stream << "-----------------------------" << std::endl;
    tri_log_stream.close();
    
}
*/
