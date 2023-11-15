/* Author: Daniel Duffy, University of Cambridge, <dld34@cam.ac.uk>
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, <dld34@cam.ac.uk>. All rights reserved.
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

Function to return the shortest distance between a point and a line segment. The
point and the endpoints of the line segment are the function arguments, given as
Eigen position vectors. This is a 3D version but exactly the same algorithm
should work for 2D. In fact Eigen probably has some fancy way of letting you
do this with some kind of matrix base class.

This approach is based on Grumdrig's answer on Stackoverflow 849211.
Another good site for this kind of thing is:
http://paulbourke.net/geometry/pointlineplane/
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>
#include <cmath>


double pointToLineSegmentDistance(
    const Eigen::Vector3d &point,
    const Eigen::Vector3d &segEndpoint0,
    const Eigen::Vector3d &segEndpoint1){

    // Calculate segment vector.
    const Eigen::Vector3d segVec = segEndpoint1 - segEndpoint0;

    // Calculate square of segment length.
    const double sqSegLength = segVec.dot(segVec);

    // Calculate vector from point to segEndpoint0.
    const Eigen::Vector3d segEndpoint0ToPointVec = point - segEndpoint0;

    // Handle case where segment has zero length.
    if( sqSegLength == 0.0 ) return segEndpoint0ToPointVec.norm();

    /* Parametrise the line containing the segement as r = r0 + t (r1-r0). Then
    t is 0 and 1 at segEndpoint0 and segEndpoint1 respectively. The closest
    point on the (infinite) line to our point must have (r-point)
    perpendicular to the line vector. This position on the line has the
    following t value:*/
    double t_ClosestPointOnInfiniteLine = segEndpoint0ToPointVec.dot(segVec) / sqSegLength;
    /* The position vector of this position on the line would be:
    segEndpoint0 + t_ClosestPointOnInfiniteLine * (segVec);
    */

    /* If the above t value is in the range (0,1), (corresponds to a point on
    the line segment), we return the corresponding distance. Otherwise, the
    shortest point-segment distance is that to the nearest segEndpoint, so we
    return that instead.*/
    if( t_ClosestPointOnInfiniteLine <= 0.0 ){
        return segEndpoint0ToPointVec.norm();
    }
    else{
        if( t_ClosestPointOnInfiniteLine >= 1.0 ){
            return (point - segEndpoint1).norm();
        }
        else{
            return sqrt(segEndpoint0ToPointVec.squaredNorm() - t_ClosestPointOnInfiniteLine * t_ClosestPointOnInfiniteLine * sqSegLength);
        }
    }
}
