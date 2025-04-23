#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Daniel Duffy dld34

Module containing the functions required to do things once a mesh (for my c++
simulation) has been created; e.g. get programmed tensors from director
patterns, plot director patterns, write a .vtk file etc.
"""

import numpy as np
from math import sin, cos, ceil, floor, copysign
import matplotlib.pyplot as plt

from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

##############################################################################
##############################################################################

# Smoothly interpolates from 0 at x=0 to 1 at x=1.
# See Inigo Quilez for alternative smoothsteps
def smoothstep(x):
    return x**2 * (3.0 - 2.0*x)

##############################################################################
##############################################################################




"""
People sometimes say that finite element stiffness matrices are 'usually banded', with no justification
for this claim. I'm confident it boils down only to this: suppose your nodes are arranged spatially like
 0  1  2  3  4  5  6  7  8  9
10 11 12 13 14 15 16 17 18 19
20 21 22 23 24 25 26 27 28 29
where the numbers are the node labels/IDs. You can see that if your domain is about n nodes wide 
(left-to-right, n=10 in this example) and m nodes tall, then spatially neighbouring nodes will have labels
differing by ar0und n or less. But the square stiffness matrix is (n*m) x (n*m). So e.g if m=n, the length of 
a r0w of the stiffness matrix is n^2, and the non-zer0 elements on that r0w will only span about n either side 
of the diagonal, so for large n you actually do have a rather banded matrix, even though you numbered your nodes
only in a very naive way! Unstructured meshes cannot look like the above in general, but similar principles 
presumably apply, at least to some extent, and people have found significant performance benefits arising from relabelling 
nodes/edges/triangles; e.g. see 
Burgess, Giles - Renumbering unstructured grids to improve the performance of codes on hierarchical memory machines
Das et al - Design and implementation of a parallel unstructured Euler solver using software primitives.
Therefore, I do the following: Apply SciPy's reverse Cuthill-Mckee implementation to the node-node adjacency matrix,
to find a new labelling of nodes where nearby node labels correspond to spatially close nodes. Then I update the 
triangulation accordingly, and then rename the triangles (just a reshuffling of the rows of triangulation) so they
are ordered by the node label of their first vertex.
"""
def optimize_mesh_tri_and_node_labels(triangulation, nodes):
    # First argument is a num_tris x 3 matrix where each element is a node id, i.e.
    # an index into the second argument, which is a num_nodes x D matrix of positions, 
    # where the dimension D is usually 2 or 3.
    
    assert(triangulation.max() == nodes.shape[0]-1), "Error: nodes and triangulation are inconsistent."

    # First we calculate the node-node adjacency matrix. That is, the square
    # matrix where the (i,j) element equals 1 if node i and node j are connected 
    # by an edge, and equals 0 otherwise. This is an adjacency matrix for an 
    # undirected graph. It is also symmetric.
    
    # The Nth element of this list is a list of the IDs of the nodes connected to N by an edge.
    # In C++ something like this would be a pretty horrible data structure if 
    # you ever wanted to iterate through it, because a std::list is "linked" --- 
    # good for inserting elements, bad for fast lookup/iteration. However, in Python
    # standard lists are not linked, and their elements are contiguous in memory,
    # so they're more like std::vector --- O(n) insertion, but fast lookup/iteration.
    neighbour_node_lists = [ [] for dummy in range(nodes.shape[0]) ] 

    for t in range(0, triangulation.shape[0]):
        for v in range(0, 3):
            for u in range(0, 3):
                neighbour_node_lists[triangulation[t,v]].append(triangulation[t,u])

    for n in range(0, nodes.shape[0]):
        # Remove node itself from its own neighbour list
        neighbour_node_lists[n][:] = (el for el in neighbour_node_lists[n] if el != n) # stackoverflow.com/a/1157174
        # Remove any duplicates in neighbour list.
        new_list = []
        [new_list.append(el) for el in neighbour_node_lists[n] if el not in new_list]
        neighbour_node_lists[n] = new_list
    
    row_idxs = []
    col_idxs = []
    vals = []
    for n in range(0, nodes.shape[0]):
        for nn in neighbour_node_lists[n]:
            row_idxs.append(n)
            col_idxs.append(nn)
            vals.append(1)
    row_idxs = np.array(row_idxs)
    col_idxs = np.array(col_idxs)
    vals = np.array(vals)

    node_adjacency_mat = coo_matrix((vals, (row_idxs, col_idxs)), shape=(nodes.shape[0], nodes.shape[0]), dtype=int)
    node_adjacency_mat = node_adjacency_mat.tocsr()
    
    
    # Then we apply the reverse Cuthill-McKee algorithm to compute a relabelling
    # of nodes that reduces the bandwidth of that adjacency matrix, which should
    # increase the tendency for spatially nearby nodes to have similar labels.
    perm = reverse_cuthill_mckee(node_adjacency_mat, symmetric_mode=True) # The Nth element of perm equals the OLD label of the node with NEW label N.
    # If you did
    # temp = node_adjacency_mat.toarray()[perm, :]
    # temp = temp[:,perm]
    # then temp would be the new node adjacency matrix, which should have much
    # smaller bandwidth than the original.
    
    
    # Now we we shuffle the nodes array so that the new
    # node labels correspond to indices into the nodes array.
    nodes[:,:] = nodes[perm, :] # The [:,;] here is needed to modify the actual original nodes matrix outside this function.
    
    
    # Got this from stackoverflow.com/a/25535723
    new_node_labels = np.empty_like(perm) # The Nth element of this equals the NEW label of the node with OLD label N.
    new_node_labels[perm] = np.arange(perm.size)
    
    
    # Now we implement the new node labels in the triangulation.
    new_triangulation = np.empty_like(triangulation)
    for t in range(0, triangulation.shape[0]):
        for v in range(0, 3):
            new_triangulation[t,v] = new_node_labels[triangulation[t,v]]
    triangulation[:,:] = new_triangulation # The [:,;] here is needed to modify the actual original nodes matrix outside this function.
    
    # Then we come up with a whole new tri numbering, such that it also has the
    # the tendency to give spatially nearby triangles similar labels. Since 
    # we've already given the node labels the analagous tendency, I take the 
    # simple approach of just sorting the triangle labels (which are just
    # indices into the triangulation array) by the node labels of their first 
    # vertices. This is probably not the optimal triangle ordering, but I'm sure
    # it's decent; it certainly seems to give essentially just as good for the 
    # spatial distribution of tri labels as we're getting for the node labels.
    triangulation[:,:] = triangulation[triangulation[:,0].argsort()] # The [:,;] here is needed to modify the actual original nodes matrix outside this function.
    
    

def discretise_Director_Pattern_Into_Bins(dirAng, range_=(0,np.pi), **kwargs):
    """
    Function that takes either:

    i) a director angle, and the desired number of 'bins' to use via a kwarg:
    'numBins'.
    An optional argument is an angular range (tuple) that is [0,pi] (the
    default) or some subset thereof. A new angle is returned, that is the input
    director angle (modulo pi) rounded to the nearest angle from the set of bin
    midpoints (found by splitting the range_ into numBins equal intervals
    ('bins') and then finding the midpoint of each interval).

    ii) a director angle, and a 1-D numpy array holding a chosen set of bin
    end points, via a kwarg: 'binEndpoints'. The rest of the functionality is
    then matches i).
    """
    # Some simple checks on a few things. These are not exhaustive. The check
    # that director angles are all real is because complex values can creep in
    # if angles were found using an inverse trig function for example.
    assert(isinstance(range_,tuple)
    and abs(range_[0]) >= 0 and abs(range_[0]) <= np.pi and abs(range_[1]) >= 0 and abs(range_[1]) <= np.pi), "Error: Something is wrong with one of the input arguments. Aborting."
    assert(np.isrealobj(dirAng) == True), "Error: Some director angles have complex values. Aborting."
    assert(not (('numBins' in kwargs) and ('binEndpoints' in kwargs))), "Error: You can specify numBins (and a range, optionally) OR an array of bin end points, but not both!"

    # Now do some modulo arithmetic to map all angles to [0, pi], the range of
    # physically distinguishable director angles.
    dirAng = dirAng % np.pi


    # Now implement functionality i).
    if ('numBins' in kwargs):
        numBins = kwargs.get('numBins')
        assert(numBins > 0), "Error: numBins must be > 0. Aborting."

        # Create some equispaced bin endpoints spanning range_.
        binEndpointsArray = np.linspace(range_[0], range_[1], numBins+1, retstep=False)

        # Now do the rounding of dirAng to the nearest bin midpoint.
        # NB a 'faster' method here would use bisect_left, but that would be
        # less readable and possibly not even faster for small numBins.

        # Two edge cases done first, unlikely to ever actually be needed.
        if dirAng < range_[0]:
            print('An edge case occured that *maybe* should be impossible - check this.')
            return (binEndpointsArray[0] + binEndpointsArray[1])/2

        elif dirAng >= range_[1]:
            print('An edge case occured that *maybe* should be impossible - check this.')
            return (binEndpointsArray[-2] + binEndpointsArray[-1])/2

        # Now the usual case.
        else:
            for i in range(0, numBins):
                if dirAng >= binEndpointsArray[i] and dirAng < binEndpointsArray[i+1]:
                    return (binEndpointsArray[i] + binEndpointsArray[i+1])/2

    # Now implement functionality ii).
    if ('binEndpoints' in kwargs):
        binEndpoints = kwargs.get('binEndpoints')
        assert(binEndpoints.ndim == 1
        and (type(binEndpoints) is np.ndarray)), "Error: binEndpoints must be a 1-D numpy array. Aborting."
        assert(np.amax(np.fabs(binEndpoints)) <= np.pi and np.amin(np.fabs(binEndpoints)) >= 0), "Error: You're trying to bin an angular range larger than pi, which makes no physical sense. Aborting."

        # First sort bin endpoints array into ascending order.
        binEndpoints.sort()


        # First bin between the provided bin end points.
        # NB a 'faster' method here would use bisect_left, but that would be
        # less readable and possibly not even faster for small numBins.
        for i in range(0, binEndpoints.shape[0]-1):
            if dirAng >= binEndpoints[i] and dirAng < binEndpoints[i+1]:
                return (binEndpoints[i] + binEndpoints[i+1])/2

        # If the bin end points don't span the full range [0, pi], take whatever is
        # left over of that interval, and combine it into one further bin.
        # I.e. take whatever is left of the range [0,pi] that hasn't already been
        # covered and use the 'extra' potential slide to do one further bin for all
        # the director angles in this remaining set of angles. If a range is left
        # over below the smallest bin end point AND above the largest, these two
        # ranges are combined here via modulo arithmetic. Any director angle edge
        # cases missed by the procedure so far should also be caught here.
        if(dirAng < binEndpoints[0] or dirAng >= binEndpoints[-1]):
            return ((binEndpoints[0] + np.pi + binEndpoints[-1])/2.0) % np.pi
        else:
            assert(False), "Error: A director angle didn't appear to belong in an bin, so something has gone wrong. Check edge cases etc."


##############################################################################
##############################################################################



def plot_Director_Pattern(dirAngs, points, fileName = None):
    # Function that takes an array of director angles, and a corresponding
    # set of points (often the centroids of triangles in a mesh), and creates a
    # plot showing a bar for each director, with its centre at the
    # corresponding point. There is also an optional file name argument - the
    # plot is saved to this filename if this argument is given.

    numAngs = dirAngs.shape[0]

    # Some simple checks on a few things. These are not exhaustive. The check
    # that director angles are all real is because complex values can creep in
    # if angles were found using an inverse trig function for example.
    assert(dirAngs.shape == (numAngs,)
    and points.shape == (numAngs,2)), "Error: Something is wrong with one of the input arguments. Aborting."
    assert(np.isrealobj(dirAngs) == True), "Error: Some director angles have complex values. Aborting."

    dirArrowXComps = np.cos(dirAngs)
    dirArrowYComps = np.sin(dirAngs)

    fig, ax = plt.subplots()
    ax.quiver(points[:,0], points[:,1], dirArrowXComps, dirArrowYComps, width=0.001, headlength=0, headaxislength=0, pivot='middle')
    ax.set_aspect('equal')
    # To turn axes off entirely:
    #ax.set_axis_off()


    if fileName is None:
        fig.show()
    else:
        assert(isinstance(fileName,str)), "Error: The file name you gave is not a string. Aborting."
        fig.savefig(fileName, format="pdf", bbox_inches = 'tight')

    return



##############################################################################
##############################################################################



def prog_Tensors_From_Top_And_Bottom_Director_Patterns(top_DirAngs, bott_DirAngs, Lambda, nu_Therm, **kwargs):
    # Function that takes the top and bottom surface director fields on a sheet,
    # some LCE material parameters, and optionally a thickness. Omitting the
    # thickness argument in calling the function sets the programmed second
    # fundamental form to zero, and this should be done in cases where there is
    # no through-thickness director variation. The function returns the
    # programmed metric and second fundamental form (secFF) components in the
    # following format: numTris x 3 arrays where each row holds the top-left,
    # then off-diagonal, then bottom right components of the programmed tensor
    # for that triangle. This is because both programmed tensors are symmetric
    # and so only have 3 independent components. Full 2x2 matrices are used
    # within the function for convenience however. The function also returns a
    # vector of numTris elements holding the values of 'tau'; the programmed
    # scalar that goes with the programmed metric.

    numTris = top_DirAngs.shape[0]

    # Some simple checks on a few things. These are not exhaustive. The check
    # that director angles are all real is because complex values can creep in
    # if angles were found using an inverse trig function for example.
    assert(top_DirAngs.shape == bott_DirAngs.shape), "Error: Something is wrong with one of the input arguments. Aborting."
    assert(np.isrealobj(top_DirAngs) == True
    and np.isrealobj(bott_DirAngs) == True), "Error: Some director angles have complex values. Aborting."

    # Find midplane director angles, assuming linear variation through
    # thickness, as is experimentally expected. Note if the top and bottom
    # directors are perpendicular to each other at any point, there is an
    # ambiguity in the midplane angle. There are two possibilities in that case
    # and this function picks one of them arbitrarily.
    midplaneDirAngs = (top_DirAngs + bott_DirAngs) / 2

    # Find the differences between the top and bottom director angles.
    dirAngDiffs = top_DirAngs - bott_DirAngs

    # Calculate a numTrisx2x2 array of rotation matrices corresponding to the
    # midplane director angles.
    rotMats = np.zeros((numTris,2, 2))
    rotMats[:,0,0] = np.cos(midplaneDirAngs)
    rotMats[:,0,1] = -np.sin(midplaneDirAngs)
    rotMats[:,1,0] = -rotMats[:,0,1]
    rotMats[:,1,1] = rotMats[:,0,0]

    ### METRIC #############
    progMetrics = np.zeros((numTris,2, 2))
    taus = np.ones(numTris) # Defauly (no twist) value is 1

    # First get principal frame values.
    # If thickness argument given, use twist nematic form, ASSUMING NU_THERM = 1/2!!!
    if ('thickness' in kwargs):
        #print('Using form for programmed metric with director twist and nu_Therm = 0.5')

        tempMats = np.zeros((numTris,2, 2))
        # Preliminary calculation.
        tempMats[:,0,0] = 2.0 * Lambda**2 / (1 + Lambda**3 + (1 - Lambda**3) * np.sinc(dirAngDiffs/np.pi))
        tempMats[:,1,1] = 2.0 * Lambda**2 / (1 + Lambda**3 - (1 - Lambda**3) * np.sinc(dirAngDiffs/np.pi))
        print("CHECK THE ABOVE LINE! SHOULD THAT LAMBDA**2 Have a nu??")

        # Now find tau.
        taus = 1.0 / np.cbrt(tempMats[:,0,0] * tempMats[:,1,1] / Lambda )

        # Finally combine these to get programmed (energetically favoured)
        # metric.
        for i in range(0, numTris):
            progMetrics[i,:,:] = taus[i] * tempMats[i,:,:]


    else:
        #print('Using no-twist form for programmed metric')
        progMetrics[:,0,0] = Lambda**2
        progMetrics[:,1,1] = Lambda**(-2*nu_Therm)


    # Now rotate into actual orientation according to midplane angle.
    for i in range(0, numTris):
        progMetrics[i,:,:] = rotMats[i,:,:] @ progMetrics[i,:,:] @ np.transpose(rotMats[i,:,:])

    # Now extract the 3 independent components of each programmed metric. These
    # will be returned by the function in a numTris x 3 array.
    progMetricComps = np.stack((progMetrics[:,0,0], progMetrics[:,0,1], progMetrics[:,1,1]), 1)


    ### SECOND FUNDAMENTAL FORM #############
    # First get principal frame values.
    # If thickness argument given, use twist nematic form, ASSUMING NU_THERM = 1/2!!!
    if ('thickness' in kwargs):
        #print('Using form for programmed secFF with director twist and nu_Therm = 0.5')

        thickness = kwargs.get('thickness')

        progSecFFs = np.zeros((numTris,2, 2))

        tempMat = np.array([[0,1],[1,0]])

        # Now apply formula for prog secFF.
        for i in range(0, numTris):
            metDetFac = progMetrics[i,0,0] * progMetrics[i,1,1] - progMetrics[i,1,0] * progMetrics[i,0,1]

            # Use Taylor expansion for small twist to avoid expression blowing up
            if abs(dirAngDiffs[i]) < 0.01:
                progSecFFs[i,:,:] = metDetFac**2.5 * (1-Lambda**3) * dirAngDiffs[i] * tempMat / (2 * Lambda**3 * thickness)

            # Use full formula for twists that are not too small
            else:
                #print(3 * metDetFac**2.5 * (1-Lambda**3) * (sin(dirAngDiffs[i]) - dirAngDiffs[i]*cos(dirAngDiffs[i])) * tempMat / (2 * Lambda**3 * thickness * dirAngDiffs[i]**2))
                progSecFFs[i,:,:] = 3 * metDetFac**2.5 * (1-Lambda**3) * (sin(dirAngDiffs[i]) - dirAngDiffs[i]*cos(dirAngDiffs[i])) * tempMat / (2 * Lambda**3 * thickness * dirAngDiffs[i]**2)

        # Now rotate into actual orientation according to midplane angle.
        for i in range(0, numTris):
            progSecFFs[i,:,:] = rotMats[i,:,:] @ progSecFFs[i,:,:] @ np.transpose(rotMats[i,:,:])

        # Now extract again the 3 independent components.
        progSecFFComps = np.stack((progSecFFs[:,0,0], progSecFFs[:,0,1], progSecFFs[:,1,1]), 1)

    # If no thickness argument supplied, programmed secFF should be zero.
    else:
        #print('Setting all programmed secFF components to zero')
        progSecFFComps = np.zeros((numTris,3))

    ###### RETURN PROGRAMMED TENSORS #######
    return progMetricComps, progSecFFComps, taus


##############################################################################
##############################################################################

def find_Edges(triangulation, nodes=None, **kwargs):
    """ This function is a more up-to-date version of boundedges(p, t) from
    distmesh utils.py. It takes a triangulation and returns a numEdges x 2
    array giving the 2 node labels corresponding to each of the edges in the
    mesh defined by the triangulation. The numTris x 3 matrix giving edge
    labels for each triangle is also returned.
    A corresponding set of 2D or 3D node positions can also optionally be
    passed; in this case a second numEdges x 2[or 3] array is also returned,
    that gives the corresponding spatial edge vectors. A final optional keyword
    argument is also present: 'just_boundary'. If this is set to False
    (default), the behaviour is just as above; if it's True, the function only
    returns arrays for the boundary edges in the mesh, and furthermore (FOR 2D
    NODE POSITIONS ONLY!!) orders each row of both arrays, so that the edges
    all point in a consistent direction going around any boundary in the mesh.
    """

    # Some simple checks on a few things. These are not exhaustive.
    assert(
        triangulation.shape[1] == 3
    and triangulation.ndim == 2
    and np.amin(triangulation)==0), "Error: Something is wrong with one of the input arguments. Aborting."

    numTris = triangulation.shape[0]

    # This first edges matrix has redundant duplicated entries. We will get rid
    # of them first
    allEdges = np.vstack((
        triangulation[:,[0,1]],
        triangulation[:,[0,2]],
        triangulation[:,[1,2]]))

    # Get rid of duplicates. uniqueEdges is just allEdges with duplicate rows
    # removed. uniqueEdgeIdxs holds the indices of these in allEdges.
    # invEdgeIdxs is the array of indices (which includes duplicates) such that
    # uniqueEdges[invEdgeIdxs] = allEdges. edgeTriCounts is the vector of counts
    # numbers of times each unique edge occured in allEdges. I.e. it's a
    # vector with an element for each row of uniqueEdges, where a 1 indicates a
    # boundary edge and a 2 indicates an internal edge.
    # NB .sort() is O(NlogN) at best, so even this first step alone means that
    # this is not going to be the fastest approach for large N; the approach in
    # my c++ code for example should be faster. Speed is not crucial here
    # though, and the approach is readable and pythonic. Note, the returned
    # uniqueEdges will be sorted vertically so that the first elements from
    # each row are in ascending order.
    allEdges.sort(1)
    uniqueEdges, uniqueEdgeIdxs, invEdgeIdxs, edgeTriCounts = np.unique(allEdges, return_index=True, return_inverse=True, return_counts=True, axis=0)

    if numTris < 5:
        print(uniqueEdges)
        print(invEdgeIdxs)
        print(uniqueEdges[invEdgeIdxs,:])
    # Now get numTris x 3 matrix where each row gives the 3 edge labels for
    # that triangle.
    triEdgeLabels = np.transpose(np.vstack((invEdgeIdxs[0:numTris],invEdgeIdxs[numTris:2*numTris],invEdgeIdxs[2*numTris:3*numTris])))

    # If you wanted to store the labels of the 1 or 2 tris sharing each edge,
    # that could be done as follows: Take triEdgeLabels and split it into its
    # 3 columns. Take the first column, and call np.arg_sort on it. Then also
    # index it with the result to get a sorted version of that column. Then
    # each element of the sorted version contains a (unique) edge label, and
    # the corresponding element of the arg_sort output vector is the label of
    # a triangle having that edge. Then do the same for the other two columns
    # of triEdgeLabels. Then think of some sensible data structure to combine
    # and store this information! Perhaps a numEdges x 2 matrix, where the
    # boundary edges have a -1 on their row as well as a single triangle label?

    # Handle case where the function will return just boundary edges.
    if kwargs.get('just_boundary') == True:

        assert(nodes is not None
        and nodes.shape[1]==2), "Error: You need to supply an array of 2D node positions to find ORIENTATED boundary edges. Aborting."

        # Find boundary edges, i.e. those that had an edgeTriCount of 1.
        boundaryUniqueEdges = uniqueEdges[edgeTriCounts==1,:]

        # otherNodes is a vector, with one element for each edge in
        # allEdges, whose contents is the label of the node in the triangle
        # that corresponds to that entry in edges, which is not either of
        # the 2 nodes in the corresponding edge.
        otherNodes = np.hstack((triangulation[:,2],triangulation[:,1],triangulation[:,0]))
        # Now retain just the entries in this corresponding to boundary
        # edges
        boundaryOtherNodes = otherNodes[uniqueEdgeIdxs[edgeTriCounts==1]]

        # Now fix edge orientations to be consistent around a boundary.
        # To do this, we look at each boundary edge, and another vector from
        # the second node of the edge to the 'other' node in the triangle, and
        # find the sign of the z component of the cross product. If it's
        # positive, the nodes in the edge are swapped.
        edgeVec = nodes[boundaryUniqueEdges[:,1]] - nodes[boundaryUniqueEdges[:,0]]
        otherVec = nodes[boundaryOtherNodes] - nodes[boundaryUniqueEdges[:,1]]
        toSwap = ( edgeVec[:,0] * otherVec[:,1] - edgeVec[:,1] * otherVec[:,0] ) > 0
        boundaryUniqueEdges[toSwap,:] = boundaryUniqueEdges[toSwap, ::-1]

        # Finally, calculate the new edge vectors after swapping, and
        # and return these as well as boundaryUniqueEdges.
        edgeVec = nodes[boundaryUniqueEdges[:,1]] - nodes[boundaryUniqueEdges[:,0]]

        return boundaryUniqueEdges, triEdgeLabels, edgeVec

    # Handle case where the function will return all unique edges, not just
    # those on the boundary.
    else:
        # We already have uniqueEdges, so we either just return that, or if
        # an array of node positions was supplied, we find the actual edge
        # vectors (no particular orientation), and then return both arrays.
        if nodes is None:

            return uniqueEdges, triEdgeLabels

        else:

            edgeVec = nodes[uniqueEdges[:,1],:] - nodes[uniqueEdges[:,0],:]

            return uniqueEdges, triEdgeLabels, edgeVec



##############################################################################
##############################################################################


def write_VTK(fileName, headerText, nodes, triangulation, abar_info, bbar_info, ref_thicknesses, ref_shear_moduli, tri_tags, constraint_indicators, node_tags):
    # Function to take the programmed tensors, and some other parameters, and
    # write these to an output file, in a (subset of) the .vtk legacy format,
    # which is the input format for my c++ simulation code, and can also be
    # read directly into ParaView.
    # nodes holds the node positions in a numNodes x 2 numpy array.
    # triangulation holds the triangulation (based on node and triangle labels
    # starting at zero) in numTris x 3 numpy array.

    # Set desired format and precision of floating point quantities.
    fp_fmt = '%.15e' #'%.5e'

    numNodes = nodes.shape[0]
    numTris = triangulation.shape[0]

    # Some simple checks on a few things. These are not exhaustive.
    assert(isinstance(fileName, str)
    and isinstance(headerText, str)
    and nodes.shape[1] == 2 or nodes.shape[1] == 3
    and triangulation.shape[1] == 3
    and constraint_indicators.shape[0] == numNodes
    and node_tags.shape[0] == numNodes), "Error: Something is wrong with one of the input arguments to write_VTK. Aborting."

    # Add 3rd column to node coordinates if necessary, to represent 3D coords.
    if nodes.shape[1] == 2:
        nodes3D = np.column_stack((nodes, np.zeros(numNodes)))
    else:
        nodes3D = nodes

    with open(fileName, 'w') as file:

        # Write preamble and node coordinates to file.
        preamble = "# vtk DataFile Version 2.0\n" + headerText + "\n" + "ASCII\n" + "DATASET POLYDATA\n" + "POINTS %d double" % numNodes
        np.savetxt(file, nodes3D, fmt=fp_fmt, delimiter=' ', header=preamble, comments='')

        # Now do triangulation with its preamble. First column of '3's is part
        # of format, indicating that each element is a triangle.
        preamble = "POLYGONS %d %d" % (numTris, 4*numTris)
        np.savetxt(file, np.column_stack((3*np.ones(numTris, dtype=int), triangulation)), fmt='%d', delimiter=' ', header=preamble, comments='')

        # Now do triangle quantities like abar_info, ref_thickness, etc, and their
        # preamble.
        num_tri_quants = 5
        file.write("CELL_DATA %d\nFIELD tri_quantities %d\n" % (numTris, num_tri_quants))

        preamble = "abar_info %d %d double" % (abar_info.shape[1], numTris)
        np.savetxt(file, abar_info, fmt=fp_fmt, delimiter=' ', header=preamble, comments='')
        preamble = "bbar_info %d %d double" % (bbar_info.shape[1], numTris)
        np.savetxt(file, bbar_info, fmt=fp_fmt, delimiter=' ', header=preamble, comments='')
        preamble = "ref_thicknesses 1 %d double" % numTris
        np.savetxt(file, ref_thicknesses, fmt=fp_fmt, delimiter=' ', header=preamble, comments='')
        preamble = "ref_shear_moduli 1 %d double" % numTris
        np.savetxt(file, ref_shear_moduli, fmt=fp_fmt, delimiter=' ', header=preamble, comments='')
        preamble = "tri_tags 1 %d int" % numTris
        np.savetxt(file, tri_tags, fmt='%d', delimiter=' ', header=preamble, comments='')
        
        # Now do node quantities like constraint indicators and other misc tags.
        num_node_quants = 2
        file.write("POINT_DATA %d\nFIELD node_quantities %d\n" % (numNodes, num_node_quants))
        preamble = "constraint_indicators 1 %d int" % numNodes
        np.savetxt(file, constraint_indicators, fmt='%d', delimiter=' ', header=preamble, comments='')
        preamble = "node_tags 1 %d int" % numNodes
        np.savetxt(file, node_tags, fmt='%d', delimiter=' ', header=preamble, comments='')


    """
    An old snippet
    # Create or open file, and write preamble.
    with open(fileName, 'w') as file:
        preamble = ["# This file is compliant with vtk DataFile Version 2.0", headerText, "ASCII", "DATASET POLYDATA", "POINTS %d double" % numNodes]
                    file.writelines("%s\n" % line for line in preamble)
    """
    return


##############################################################################
##############################################################################

def read_Between_Lines(fileName, substr1, substr2 = None, **kwargs):
    # This function returns a block of text/data taken from the file specified
    # in fileName, either as a list of lines, or as a numpy array. The first
    # block is found in the file that is between a line containing substr1 and a
    # line containing substr2 if supplied. The start and end ines are NOT
    # included in the block. Optionally, the final argument substr2 can be
    # omitted, in which case the end of the block will be the end of the file.
    # This will also happen if substr1 is found but no instance of substr2 is
    # found thereafter. All leading and trailing whitespace (including
    # newlines \n) is stripped from the lines to be returned. Additionally a int
    # kwarg skip_lines can be used to delete this many lines from the start of
    # the block before returning. A final kwarg can also be given, indicating
    # that a numpy array should be returned instead of a list of lines, where
    # each line is split into its 'words' to create each line of the numpy array.
    # The data type of the numpy array is specified by the kwarg and should be a
    # string. So for example, the function's final argument could be:
    # return_np_array_using_dtype = 'float'.
    # If returning a numpy array, for now any commas present on each line are
    # replaced by spaces, which is reasonable behaviour for most scientific
    # applications.
    # More complex versions of this function could be constructed easily with
    # the Regex module; see https://www.computerhope.com/issues/ch001721.htm
    # for inspiration. I made the matching for substr1 and substr2 case
    # insensitive, but that could easily be changed.

    # Some simple checks on a few things. These are not exhaustive.
    assert(isinstance(fileName, str)
    and isinstance(substr1, str)
    and (substr2 == None or isinstance(substr2, str))), "Error: Something is wrong with one of the input arguments to read_Between_Lines. Aborting."

    linesToReturn = []

    # Make the two substrings case insensitive.
    substr1 = substr1.casefold()
    if substr2 != None: substr2 = substr2.casefold()

    with open(fileName, 'r') as file:

        inDesiredBlock = False

        for line in file:

            # Test whether we have reached the start or end of the desired
            # block, and set inDesiredBlock accordingly.
            # if line.casefold().find(substr1) != -1: inDesiredBlock = True # older version
            if substr1 in line.casefold():
                inDesiredBlock = True
                continue
            if substr2 != None:
                if substr2 in line.casefold():
                    inDesiredBlock = False
                    # Terminate once reached end of desired block.
                    break

            # Copy current line if in desired block, stripping leading and
            # trailing whitespace.
            if inDesiredBlock == True:  linesToReturn.append(line.strip())

    # Remove some lines from start of linesToReturn if skip_lines is specified,
    if kwargs.get('skip_lines'):
        skip_lines = kwargs.get('skip_lines')
        assert(isinstance(skip_lines, int)
        and (not skip_lines < 1)), "Error: skip_lines must be an int > 0. Aborting."

        del linesToReturn[0:skip_lines]

    # Return np.array if specified.
    if kwargs.get('return_np_array_using_dtype'):

        dataType = kwargs.get('return_np_array_using_dtype')

        assert(isinstance(dataType, str)), "Error: To return a numpy array you must use e.g. return_np_array_using_dtype = 'float' where 'float' can be any string specifying a numpy data type. Aborting."

        # Split each line, following Stackoverflow 4081217, also replacing any
        # comma delimeters with spaces.
        for idx, element in enumerate(linesToReturn):
            linesToReturn[idx] = element.replace(',',' ').split()

        # Return numpy array of the linesToReturn we've ended up with, removing
        # any superfluous dimensions for simplicity and intuitiveness.
        return np.squeeze(np.array(linesToReturn, dtype=dataType))

    # Otherwise just return list of lines.
    else:
        return linesToReturn

#############################################################################
#############################################################################


def get_map2DPointToSurface_Function(refNodes, defNodes, triangulation):
    # Function that returns a function to map 2D points to 3D space based on
    # the piecewise linear mapping defined by the reference and deformed states
    # of a triangulated mesh.

    numNodes = refNodes.shape[0]
    numTris = triangulation.shape[0]

    # There may or may not be a slightly faster though less intuitive way
    # to do this using barycentric coords - see various answers on
    # Stackoverflow 2049582, where I also found the method used here.
    # It works as follows: pick any vertex and draw the two sides outwards
    # from it. For the point to be in the triangle, the the corresponding
    # 'cross-products' must have opposite sign. Then look at the sign of
    # the final side's cross product with a vector from a different vertex
    # to the point to finish the determination.
    # I think in rare cases where a point is almost exactly on an edge,
    # this method can fail (see
    # http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html)
    # and you'll have to do more work in that case. For now I just check
    # that there is no such problem i.e. all nodes are assigned to at least
    # one triangle. Something you could do to detect the edge case is catch
    # when an edge think a point in on its left (say) when going along the
    # edge in either direction. Then you could pick some convention etc to
    # handle it.
    def isPointInTriangle(point, vert0, vert1, vert2):

        vert0ToPoint_x = point[0] - vert0[0]
        vert0ToPoint_y = point[1] - vert0[1]

        check1 = (vert1[0]-vert0[0]) * vert0ToPoint_y - (vert1[1]-vert0[1]) * vert0ToPoint_x > 0

        if ( (vert2[0]-vert0[0]) * vert0ToPoint_y - (vert2[1]-vert0[1]) * vert0ToPoint_x > 0 ) == check1:

            return False

        if ( (vert2[0]-vert1[0]) * (point[1] - vert1[1]) - (vert2[1]-vert1[1]) * (point[0] - vert1[0]) > 0 ) != check1:

            return False

        return True

    class Tri:
       def  __init__(self, vertLabels):
           self.vertLabels = vertLabels
           self.defVerts = -9876*np.ones((3,3), dtype = 'float')
           self.refVerts = -4321*np.ones((2,2), dtype = 'float')
           self.defSides = -999*np.ones((3,2), dtype = 'float')
           self.refSides = -9999*np.ones((2,2), dtype = 'float')
           self.defGrad = -99999*np.ones((3,2), dtype = 'float')
           self.refCentroid = -12345.6*np.ones(2)
           return

       def set_defVerts(self, defVerts):
           self.defVerts = defVerts
           return

       def set_refVerts(self, refVerts):
           self.refVerts = refVerts
           return

       def set_refCentroid(self, refCentroid):
           self.refCentroid = refCentroid
           return

       def reorderVerts_And_Calc_Sides_And_DefGrad(self):
           self.refSides[:,0] = self.refVerts[1,:] - self.refVerts[0,:]
           self.refSides[:,1] = self.refVerts[2,:] - self.refVerts[0,:]

           # Fix orientation to match c++ code
           if self.refSides[0,0] * self.refSides[1,1] - self.refSides[0,1] * self.refSides[1,0] < 0:

               temp = self.vertLabels[2]
               self.vertLabels[2] = self.vertLabels[1]
               self.vertLabels[1] = temp
               self.refVerts[0,:]

               temp2 = np.copy(self.refVerts[2,:])
               self.refVerts[2,:] = self.refVerts[1,:]
               self.refVerts[1,:] = temp2

               temp3 = np.copy(self.defVerts[2,:])
               self.defVerts[2,:] = self.defVerts[1,:]
               self.defVerts[1,:] = temp3

               self.refSides[:,0] = self.refVerts[1,:] - self.refVerts[0,:]
               self.refSides[:,1] = self.refVerts[2,:] - self.refVerts[0,:]


           self.defSides[:,0] = self.defVerts[1,:] - self.defVerts[0,:]
           self.defSides[:,1] = self.defVerts[2,:] - self.defVerts[0,:]

           self.defGrad = self.defSides @ np.linalg.inv(self.refSides)
           return

       def mapRefPointToDeformedState(self, refPoint):
          return self.defVerts[0,:] + self.defGrad @ (refPoint - self.refVerts[0,:]) # Could use any vertex


    triangles = np.empty(numTris, dtype=object)
    for t in range(0, numTris):
        triangles[t] = Tri(triangulation[t,:])
        triangles[t].set_defVerts(defNodes[triangulation[t,:]])
        triangles[t].set_refVerts(refNodes[triangulation[t,:]])
        triangles[t].set_refCentroid((refNodes[triangulation[t,0],:] + refNodes[triangulation[t,1],:] + refNodes[triangulation[t,2],:]) / 3)
        triangles[t].reorderVerts_And_Calc_Sides_And_DefGrad()


    minRefEdgeLength = np.linalg.norm(triangles[0].refSides[:,0])
    for t in range(0, numTris):
        edgeLength = np.linalg.norm(triangles[t].refSides[:,0])
        if minRefEdgeLength > edgeLength:
            minRefEdgeLength = edgeLength
        edgeLength = np.linalg.norm(triangles[t].refSides[:,1])
        if minRefEdgeLength > edgeLength:
            minRefEdgeLength = edgeLength



    gridSquareSize = minRefEdgeLength # This I'm sure is only roughly optimal at best, and could be tweaked.
    # The 2*minRefEdgeLength's are largely cautionary, and may be unnecessary overkill.
    grid_xmin = np.amin(refNodes[:,0]) - 2 * minRefEdgeLength
    grid_xmax_approx = np.amax(refNodes[:,0]) + 2 * minRefEdgeLength
    gridPoint_xvals = np.linspace( grid_xmin, grid_xmin + ceil((grid_xmax_approx-grid_xmin)/gridSquareSize)*gridSquareSize, ceil((grid_xmax_approx-grid_xmin)/gridSquareSize)+1 )
    grid_xmax = np.amax(gridPoint_xvals)

    grid_ymin = np.amin(refNodes[:,1]) - 2 * minRefEdgeLength
    grid_ymax_approx = np.amax(refNodes[:,1]) + 2 * minRefEdgeLength
    gridPoint_yvals = np.linspace( grid_ymin, grid_ymin + ceil((grid_ymax_approx-grid_ymin)/gridSquareSize)*gridSquareSize, ceil((grid_ymax_approx-grid_ymin)/gridSquareSize)+1 )
    grid_ymax = np.amax(gridPoint_yvals)


    xx, yy = np.meshgrid(gridPoint_xvals, gridPoint_yvals, sparse=False, indexing='ij') # Matrices of x and y grid point positions respectively



    class Square:
       def  __init__(self, vertPoints):
           self.Idxs = []
           self.vertPoints = vertPoints
           sides = np.zeros((2,4))
           for s in range(0,4):
               sides[:,s] = vertPoints[(s+1)%4,:] - vertPoints[s,:]
           self.sides = sides
           self.triIntersections = []
           return

       def set_Idxs(self, Idxs):
           self.Idxs = Idxs
           return

       def append_TriIntersection(self, triLabel):
           if (not (triLabel in self.triIntersections)): # if statement is precautionary, should never be necessary in this code.
               self.triIntersections.append(triLabel)
           return

    squares = np.empty((xx.shape[0]-1,xx.shape[1]-1), dtype=object)
    for i in range(0, squares.shape[0]):
        for j in range(0, squares.shape[1]):
            squares[i,j] = Square(np.array([[xx[i,j], yy[i,j]], [xx[i+1,j], yy[i+1,j]], [xx[i+1,j+1], yy[i+1,j+1]], [xx[i,j+1], yy[i,j+1]]]))
            squares[i,j].set_Idxs([i,j])
    numSquares = squares.size
    print('Number of squares in grid for get_map2DPointToSurface_Function is ' +str(numSquares))

    def SquareIdxsOfPoint(point):
        # Returns the i,j indexes in the squares array corresponding to the square
        # that a given x,y point is in.
        # assert(point[0] >= grid_xmin and point[1] >= grid_ymin and point[0] <= grid_xmax and point[1] <= grid_ymax), 'The point you provide to SquareIdxsOfPoint() must lie within the grid of squares!'
        if( point[0] >= grid_xmin and point[0] <= grid_xmax  and point[1] >= grid_ymin and point[1] <= grid_ymax ): # if point lies within the grid of squares
            return (floor((point[0]-grid_xmin) / gridSquareSize),  floor((point[1]-grid_ymin) / gridSquareSize))
        else:
            return False


    def convexPolyIntersectionCheck(points1, edges1, points2, edges2):
        # Function that returns bool indicating whether two convex 2D polygons
        # intersect/overlap, based on separating axis theorem. Each row of the
        # points1,2 inputs should hold point coords, going round anticlockwise, and
        # the nth column of edges1,2 should be a side vector from the nth point to
        # the n+1th (so also anticlockwise).
        # See Ericson 'Real-Time Collision Detection'.
        # Some optimisation of this might be possible, see e.g. http://web.archive.org/web/20141127210836/http://content.gpwiki.org/index.php/Polygon_Collision

        # I have not thought about edge cases here where the polygons
        # exactly share a point or nearly do and floating point causes problems etc,
        # as it really shouldn't matter for the current application.
        pointSepVecs = np.empty((points1.shape[0], points2.shape[0], 2))
        for p1 in range(0, points1.shape[0]):
            for p2 in range(0, points2.shape[0]):
                    pointSepVecs[p1,p2,:] = points2[p2,:] - points1[p1,:]


        # Check if any of polygon 1's edges give a separating axis
        for e in range(0, edges1.shape[1]):

            abandonThisEdge = False

            for p in range(0, points2.shape[0]):

                # For an edge to give a separating axis, all the other polygon's
                # points must lie to the right of the edge as you going along it
                # according to our anticlockwise convention.
                # If this edge does not give a separating axis, move on.
                if int(copysign(1, pointSepVecs[e,p,0] * edges1[1,e] - pointSepVecs[e,p,1] * edges1[0,e])) == -1:

                    abandonThisEdge = True
                    break

            if abandonThisEdge == True:
                continue
            # If a separating axis is found, the polygons do not intersect.
            return False

         # Check if any of polygon 2's edges give a separating axis
        for e in range(0, edges2.shape[1]):

            abandonThisEdge = False

            for p in range(0, points1.shape[0]):

                # If this edge does not give a separating axis, move on.
                if int(copysign(1, pointSepVecs[p,e,0] * edges2[1,e] - pointSepVecs[p,e,1] * edges2[0,e])) == 1:

                    abandonThisEdge = True
                    break

            if abandonThisEdge == True:
                continue

            return False

        # If no separating axis found, the polygons intersect
        return True

    # Let the squares know which triangles intersect them.
    for t in range(0, numTris):
        min_x_refVert = triangles[t].refVerts[ np.argmin(triangles[t].refVerts[:,0]), :]
        max_x_refVert = triangles[t].refVerts[ np.argmax(triangles[t].refVerts[:,0]), :]
        min_y_refVert = triangles[t].refVerts[ np.argmin(triangles[t].refVerts[:,1]), :]
        max_y_refVert = triangles[t].refVerts[ np.argmax(triangles[t].refVerts[:,1]), :]

        allThreeRefSides = np.column_stack((triangles[t].refSides[:,0], triangles[t].refSides[:,1] - triangles[t].refSides[:,0], -triangles[t].refSides[:,1]))

        # We test all squares in the grid-aligned rectangle defined by the tri verts
        idxsToCheck_x =  range(SquareIdxsOfPoint(min_x_refVert)[0], SquareIdxsOfPoint(max_x_refVert)[0]+1)
        idxsToCheck_y =  range(SquareIdxsOfPoint(min_y_refVert)[1], SquareIdxsOfPoint(max_y_refVert)[1]+1)

        for i in idxsToCheck_x:
            for j in idxsToCheck_y:

                if convexPolyIntersectionCheck(triangles[t].refVerts, allThreeRefSides, squares[i,j].vertPoints, squares[i,j].sides) == True:
                    squares[i,j].append_TriIntersection(t)



    def map2DPointToSurface(point):
        squarePointIsInIfAny = SquareIdxsOfPoint(point)
        # If point is in a square...
        if squarePointIsInIfAny != False:

            # If point's square intersects any triangles...
            possibleTris = squares[squarePointIsInIfAny].triIntersections

            for t in possibleTris:
                # If point is in THIS triangle, record that and its plane->plane mapping will be used to map the point.
                if isPointInTriangle(point,  triangles[t].refVerts[0,:], triangles[t].refVerts[1,:], triangles[t].refVerts[2,:]) == True:

                    return triangles[t].mapRefPointToDeformedState(point)

        # Otherwise (i.e. the square had no triangle
        # intersections, or it did but the point wasn't in any of those triangles),
        # we instead just find the closest square to the point that does have
        # intersections, and use the triangle it intersects whose centroid is
        # closest to the point. You could be more careful and really find e.g. the
        # closest point in the set of triangles to the given point, and use the
        # corresponding triangle. However the approach here should be fine for
        # our application. For now we exclude cases where the point fell outside
        # the square grid entirely, but that should be fixed in future, by imagining
        # the grid did extend out to where the point is. The expanding-ring-of-squares
        # search is also done a bit wastefully- in particular when it reaches the edge
        # of the grid in one direction it could just stop rather than checking new
        # 'possibleSquare's in that direction for potential adding to nextSquaresToSearch.
        if squarePointIsInIfAny == False:
                assert(False), 'Error - for now we need all points you want to map to be within the square grid, so make the grid bigger if need be.'
        else:
            possibleTris = []
            squaresToSearch = [squarePointIsInIfAny]
            startSquare = squarePointIsInIfAny
            ringIdxRad = 0 # Search based on a square 'ring' of squares, with the oouter squares this many indices from the central start square

            while True:

                for sq in squaresToSearch:
                    possibleTris.extend(squares[sq].triIntersections)

                if len(possibleTris) != 0:
                    break

                nextSquaresToSearch = []
                for sq in squaresToSearch:

                    if sq[0] == startSquare[0] - ringIdxRad:
                        nextSquaresToSearch.append((sq[0]-1, sq[1]))

                    elif sq[0] == startSquare[0] + ringIdxRad:
                        nextSquaresToSearch.append((sq[0]+1, sq[1]))

                    elif sq[1] == startSquare[1] - ringIdxRad:
                        nextSquaresToSearch.append((sq[0], sq[1]-1))

                    elif sq[1] == startSquare[1] + ringIdxRad:
                        nextSquaresToSearch.append((sq[0], sq[1]+1))

                # Corners of ring add an extra square
                nextSquaresToSearch.append((startSquare[0] + ringIdxRad + 1, startSquare[1] + ringIdxRad + 1))
                nextSquaresToSearch.append((startSquare[0] - ringIdxRad - 1, startSquare[1] + ringIdxRad + 1))
                nextSquaresToSearch.append((startSquare[0] - ringIdxRad - 1, startSquare[1] - ringIdxRad - 1))
                nextSquaresToSearch.append((startSquare[0] + ringIdxRad + 1, startSquare[1] - ringIdxRad - 1))

                # Eliminate any nextSquaresToSearch that don't exist as they're out of the grid.
                realNextSquaresToSearch = [nsq for nsq in nextSquaresToSearch if (nsq[0] >= 0 and nsq[0] < squares.shape[0] and nsq[1] >= 0 and nsq[1] < squares.shape[1])]

                squaresToSearch = realNextSquaresToSearch

                ringIdxRad += 1

            # Now we've found at least one possible tri, find the one with the closest
            # centroid (probably not the best approach but simple and should work
            # fine for now).
            centroidDists = np.array([np.linalg.norm(triangles[t].refCentroid - point) for t in possibleTris])
            triToUse = possibleTris[np.argmin(centroidDists)]

            return triangles[triToUse].mapRefPointToDeformedState(point)

        assert(False), 'Error - the map2DPointToSurface function has gone wrong somehow and a node has slipped through and not been assigned a triangle.'
        return

    # Now return the mapping function we've generated.
    return map2DPointToSurface



#############################################################################
#############################################################################


def get_piecewise_linear_interpolation_function_for_tri_mesh(nodes, nodal_field_vals, triangulation):
    # You feed this function a set of nodes and a triangulation that define a triangle
    # mesh, as well as a set of values f_i representing the values of some quantity f
    # at those nodes. This function returns a function that linearly interpolates
    # f between the nodes of the triangle mesh. So it's a piecewise-linear approximation
    # to the underlying field f(x,y). So using this function usually looks like this:
    # linterp = get_piecewise_linear_interpolation_function_for_tri_mesh(nodes, nodal_temperatures, triangulation)
    # for point in list_of_xy_points_that_are_not_nodes_of_the_mesh:
     #   temperature[j] = linterp(point)
    # The points at which you want to evaluate the linear interpolation should
    # be either within the boundaries of the mesh or very close to the boundary.
    # (This function is almost a copy-paste of get_map2DPointToSurface_Function.)
    # I expect this function could be made a lot faster, possibly by removing some
    # of the class structure.


    num_nodes = nodes.shape[0]
    num_tris = triangulation.shape[0]

    # There may or may not be a slightly faster though less intuitive way
    # to do this using barycentric coords - see various answers on
    # Stackoverflow 2049582, where I also found the method used here.
    # It works as follows: pick any vertex and draw the two sides outwards
    # from it. For the point to be in the triangle, the the corresponding
    # 'cross-products' must have opposite sign. Then look at the sign of
    # the final side's cross product with a vector from a different vertex
    # to the point to finish the determination.
    # I think in rare cases where a point is almost exactly on an edge,
    # this method can fail (see
    # http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html)
    # and you'll have to do more work in that case. For now I just check
    # that there is no such problem i.e. all nodes are assigned to at least
    # one triangle. Something you could do to detect the edge case is catch
    # when an edge think a point in on its left (say) when going along the
    # edge in either direction. Then you could pick some convention etc to
    # handle it.
    def isPointInTriangle(point, vert0, vert1, vert2):

        vert0ToPoint_x = point[0] - vert0[0]
        vert0ToPoint_y = point[1] - vert0[1]

        check1 = (vert1[0]-vert0[0]) * vert0ToPoint_y - (vert1[1]-vert0[1]) * vert0ToPoint_x > 0

        if ( (vert2[0]-vert0[0]) * vert0ToPoint_y - (vert2[1]-vert0[1]) * vert0ToPoint_x > 0 ) == check1:

            return False

        if ( (vert2[0]-vert1[0]) * (point[1] - vert1[1]) - (vert2[1]-vert1[1]) * (point[0] - vert1[0]) > 0 ) != check1:

            return False

        return True

    class Tri:
       def  __init__(self, vertex_labels):
           self.vertex_labels = vertex_labels
           self.vertex_positions = -4321*np.ones((3,2), dtype = 'float')
           self.edges = 4321*np.ones((3,2), dtype = 'float')
           self.vertex_field_vals = -9876*np.ones(3, dtype = 'float')
           self.centroid = -12345*np.ones(2)
           self.coeffs_of_linear_interpolation = -12345*np.ones(3)
           return

       def set_vertex_positions(self, vertex_positions):
           self.vertex_positions = vertex_positions
           return
       
       # Reorder the storing of vertex labels and positions so the increasing
       # index from 0 to 2 corresponds to traversing the tri anticlockwise.
       def reorder_vertices(self):
           
           vec_0_to_1 = self.vertex_positions[1,:] -self.vertex_positions[0,:]
           vec_0_to_2 = self.vertex_positions[2,:] -self.vertex_positions[0,:]
           
           if vec_0_to_1[0]*vec_0_to_2[1] - vec_0_to_1[1]*vec_0_to_2[0] < 0:
               
               temp = self.vertex_labels[1]
               self.vertex_labels[1] = self.vertex_labels[2]
               self.vertex_labels[2] = temp
               
               temp2 = self.vertex_positions[1,:]
               self.vertex_positions[1,:] = self.vertex_positions[2,:]
               self.vertex_positions[2,:] = temp2
           

       def calc_edges(self):
           # The first edge is the vector pointing from vertex 0 to vertex 1. 
           # The second points from vertex 1 to vertex 2. The third points from
           # vertex 2 to vertex 0. We calculate the edges after putting the 
           # vertices in anticlockwise order, so the calculated edges will all
           # also point anticlockwise. Each ROW of the edges matrix is an edge vector.
           self.edges[0,:] = self.vertex_positions[1,:] - self.vertex_positions[0,:]
           self.edges[1,:] = self.vertex_positions[2,:] - self.vertex_positions[1,:]
           self.edges[2,:] = self.vertex_positions[0,:] - self.vertex_positions[2,:]
       
       def set_vertex_field_vals(self, vertex_field_vals):
           self.vertex_field_vals = vertex_field_vals
           return

       def set_centroid(self, centroid):
           self.centroid = centroid
           return

       def calc_coeffs_of_linear_interpolation(self):
           # Within the triangle, the linear interpolated field will be
           # f = a + b * (x-x_centroid) + c * (y-y_centroid). To find a, b, and c,
           # we require that this linear polynomial interpolates the 3 vertex values of f.
           # The three requirements are written as mat_to_invert @ [a, b, c] = vertex_field_vals.
           mat_to_invert = np.ones((3,3), dtype = float)
           mat_to_invert[:,1] = self.vertex_positions[:,0] - self.centroid[0]
           mat_to_invert[:,2] = self.vertex_positions[:,1] - self.centroid[1]
           self.coeffs_of_linear_interpolation = np.linalg.inv(mat_to_invert) @ self.vertex_field_vals
           return

       def calc_interpolated_value(self, point):
           return self.coeffs_of_linear_interpolation[0] + self.coeffs_of_linear_interpolation[1] * (point[0]-self.centroid[0]) + self.coeffs_of_linear_interpolation[2] * (point[1]-self.centroid[1]) 


    triangles = np.empty(num_tris, dtype=object)
    for t in range(0, num_tris):
        triangles[t] = Tri(triangulation[t,:])
        triangles[t].set_vertex_positions(nodes[triangulation[t,:],:])
        triangles[t].reorder_vertices()
        triangles[t].calc_edges()
        triangles[t].set_vertex_field_vals(nodal_field_vals[triangles[t].vertex_labels])
        triangles[t].set_centroid( np.sum(triangles[t].vertex_positions, axis = 0) / 3)
        triangles[t].calc_coeffs_of_linear_interpolation()


    min_edge_length = np.linalg.norm(triangles[0].edges[0,:])
    for t in range(0, num_tris):
        length = np.linalg.norm(triangles[t].edges[0,:])
        if min_edge_length > length:
            min_edge_length = length
        length = np.linalg.norm(triangles[t].edges[1,:])
        if min_edge_length > length:
            min_edge_length = length
        length = np.linalg.norm(triangles[t].edges[2,:])
        if min_edge_length > length:
            min_edge_length = length



    gridSquareSize = 1.1 * min_edge_length # This I'm sure is only roughly optimal at best, and could be tweaked.
    # The 2*min_edge_length's are largely cautionary, and may be unnecessary overkill.
    grid_xmin = np.amin(nodes[:,0]) - 2 * min_edge_length
    grid_xmax_approx = np.amax(nodes[:,0]) + 2 * min_edge_length
    gridPoint_xvals = np.linspace( grid_xmin, grid_xmin + ceil((grid_xmax_approx-grid_xmin)/gridSquareSize)*gridSquareSize, ceil((grid_xmax_approx-grid_xmin)/gridSquareSize)+1 )
    grid_xmax = np.amax(gridPoint_xvals)

    grid_ymin = np.amin(nodes[:,1]) - 2 * min_edge_length
    grid_ymax_approx = np.amax(nodes[:,1]) + 2 * min_edge_length
    gridPoint_yvals = np.linspace( grid_ymin, grid_ymin + ceil((grid_ymax_approx-grid_ymin)/gridSquareSize)*gridSquareSize, ceil((grid_ymax_approx-grid_ymin)/gridSquareSize)+1 )
    grid_ymax = np.amax(gridPoint_yvals)


    xx, yy = np.meshgrid(gridPoint_xvals, gridPoint_yvals, sparse=False, indexing='ij') # Matrices of x and y grid point positions respectively



    class Square:
       def  __init__(self, vertPoints):
           self.Idxs = []
           self.vertPoints = vertPoints
           sides = np.zeros((2,4))
           for s in range(0,4):
               sides[:,s] = vertPoints[(s+1)%4,:] - vertPoints[s,:]
           self.sides = sides
           self.triIntersections = []
           return

       def set_Idxs(self, Idxs):
           self.Idxs = Idxs
           return

       def append_TriIntersection(self, triLabel):
           if (not (triLabel in self.triIntersections)): # if statement is precautionary, should never be necessary in this code.
               self.triIntersections.append(triLabel)
           return

    squares = np.empty((xx.shape[0]-1,xx.shape[1]-1), dtype=object)
    for i in range(0, squares.shape[0]):
        for j in range(0, squares.shape[1]):
            squares[i,j] = Square(np.array([[xx[i,j], yy[i,j]], [xx[i+1,j], yy[i+1,j]], [xx[i+1,j+1], yy[i+1,j+1]], [xx[i,j+1], yy[i,j+1]]]))
            squares[i,j].set_Idxs([i,j])
    numSquares = squares.size
    print('Number of squares in grid for get_piecewise_linear_interpolation_function_for_tri_mesh is ' +str(numSquares))

    def SquareIdxsOfPoint(point):
        # Returns the i,j indexes in the squares array corresponding to the square
        # that a given x,y point is in.
        # assert(point[0] >= grid_xmin and point[1] >= grid_ymin and point[0] <= grid_xmax and point[1] <= grid_ymax), 'The point you provide to SquareIdxsOfPoint() must lie within the grid of squares!'
        if( point[0] >= grid_xmin and point[0] <= grid_xmax  and point[1] >= grid_ymin and point[1] <= grid_ymax ): # if point lies within the grid of squares
            return (floor((point[0]-grid_xmin) / gridSquareSize),  floor((point[1]-grid_ymin) / gridSquareSize))
        else:
            return False # point does not lie within the grid of squares


    def convexPolyIntersectionCheck(points1, edges1, points2, edges2):
        # Function that returns bool indicating whether two convex 2D polygons
        # intersect/overlap, based on separating axis theorem. Each ROW of the
        # points1,2 inputs should hold point coords, going round anticlockwise, and
        # the nth COLUMN of edges1,2 should be a side vector from the nth point to
        # the n+1th (so also anticlockwise).
        # See Ericson 'Real-Time Collision Detection'.
        # Some optimisation of this might be possible, see e.g. http://web.archive.org/web/20141127210836/http://content.gpwiki.org/index.php/Polygon_Collision

        # I have not thought about edge cases here where the polygons
        # exactly share a point or nearly do and floating point causes problems etc,
        # as it really shouldn't matter for the current application.
        pointSepVecs = np.empty((points1.shape[0], points2.shape[0], 2))
        for p1 in range(0, points1.shape[0]):
            for p2 in range(0, points2.shape[0]):
                    pointSepVecs[p1,p2,:] = points2[p2,:] - points1[p1,:]


        # Check if any of polygon 1's edges give a separating axis
        for e in range(0, edges1.shape[1]):

            abandonThisEdge = False

            for p in range(0, points2.shape[0]):

                # For an edge to give a separating axis, all the other polygon's
                # points must lie to the right of the edge as you going along it
                # according to our anticlockwise convention.
                # If this edge does not give a separating axis, move on.
                if int(copysign(1, pointSepVecs[e,p,0] * edges1[1,e] - pointSepVecs[e,p,1] * edges1[0,e])) == -1:

                    abandonThisEdge = True
                    break

            if abandonThisEdge == True:
                continue
            # If a separating axis is found, the polygons do not intersect.
            return False

         # Check if any of polygon 2's edges give a separating axis
        for e in range(0, edges2.shape[1]):

            abandonThisEdge = False

            for p in range(0, points1.shape[0]):

                # If this edge does not give a separating axis, move on.
                if int(copysign(1, pointSepVecs[p,e,0] * edges2[1,e] - pointSepVecs[p,e,1] * edges2[0,e])) == 1:

                    abandonThisEdge = True
                    break

            if abandonThisEdge == True:
                continue

            return False

        # If no separating axis found, the polygons intersect
        return True

    # Let the squares know which triangles intersect them.
    for t in range(0, num_tris):
        min_x_vert = triangles[t].vertex_positions[ np.argmin(triangles[t].vertex_positions[:,0]), :]
        max_x_vert = triangles[t].vertex_positions[ np.argmax(triangles[t].vertex_positions[:,0]), :]
        min_y_vert = triangles[t].vertex_positions[ np.argmin(triangles[t].vertex_positions[:,1]), :]
        max_y_vert = triangles[t].vertex_positions[ np.argmax(triangles[t].vertex_positions[:,1]), :]

        # We test all squares in the grid-aligned rectangle defined by the tri verts
        idxsToCheck_x =  range(SquareIdxsOfPoint(min_x_vert)[0], SquareIdxsOfPoint(max_x_vert)[0]+1)
        idxsToCheck_y =  range(SquareIdxsOfPoint(min_y_vert)[1], SquareIdxsOfPoint(max_y_vert)[1]+1)

        for i in idxsToCheck_x:
            for j in idxsToCheck_y:

                if convexPolyIntersectionCheck(triangles[t].vertex_positions, triangles[t].edges.transpose(), squares[i,j].vertPoints, squares[i,j].sides) == True:
                    squares[i,j].append_TriIntersection(t)



    def interpolating_function(point):
        squarePointIsInIfAny = SquareIdxsOfPoint(point)
        # If point is in a square...
        if squarePointIsInIfAny != False:

            # If point's square intersects any triangles...
            possibleTris = squares[squarePointIsInIfAny].triIntersections

            for t in possibleTris:
                # If point is in THIS triangle, record that and its plane->plane mapping will be used to map the point.
                if isPointInTriangle(point,  triangles[t].vertex_positions[0,:], triangles[t].vertex_positions[1,:], triangles[t].vertex_positions[2,:]) == True:

                    return triangles[t].calc_interpolated_value(point)

        # Otherwise (i.e. the square had no triangle
        # intersections, or it did but the point wasn't in any of those triangles),
        # we instead just find the closest square to the point that does have
        # intersections, and use the triangle it intersects whose centroid is
        # closest to the point. You could be more careful and really find e.g. the
        # closest point in the set of triangles to the given point, and use the
        # corresponding triangle. However the approach here should be fine for
        # our application. For now we exclude cases where the point fell outside
        # the square grid entirely, but that should be fixed in future, by imagining
        # the grid did extend out to where the point is. The expanding-ring-of-squares
        # search is also done a bit wastefully- in particular when it reaches the edge
        # of the grid in one direction it could just stop rather than checking new
        # 'possibleSquare's in that direction for potential adding to nextSquaresToSearch.
        if squarePointIsInIfAny == False:
                assert(False), 'Error: For now we need all points you want to map to be within the square grid, so make the grid bigger if need be.'
        else:
            possibleTris = []
            squaresToSearch = [squarePointIsInIfAny]
            startSquare = squarePointIsInIfAny
            ringIdxRad = 0 # Search based on a square 'ring' of squares, with the outer squares this many indices from the central start square

            while True:

                for sq in squaresToSearch:
                    possibleTris.extend(squares[sq].triIntersections)

                if len(possibleTris) != 0:
                    break

                nextSquaresToSearch = []
                for sq in squaresToSearch:

                    if sq[0] == startSquare[0] - ringIdxRad:
                        nextSquaresToSearch.append((sq[0]-1, sq[1]))

                    elif sq[0] == startSquare[0] + ringIdxRad:
                        nextSquaresToSearch.append((sq[0]+1, sq[1]))

                    elif sq[1] == startSquare[1] - ringIdxRad:
                        nextSquaresToSearch.append((sq[0], sq[1]-1))

                    elif sq[1] == startSquare[1] + ringIdxRad:
                        nextSquaresToSearch.append((sq[0], sq[1]+1))

                # Corners of ring add an extra square
                nextSquaresToSearch.append((startSquare[0] + ringIdxRad + 1, startSquare[1] + ringIdxRad + 1))
                nextSquaresToSearch.append((startSquare[0] - ringIdxRad - 1, startSquare[1] + ringIdxRad + 1))
                nextSquaresToSearch.append((startSquare[0] - ringIdxRad - 1, startSquare[1] - ringIdxRad - 1))
                nextSquaresToSearch.append((startSquare[0] + ringIdxRad + 1, startSquare[1] - ringIdxRad - 1))

                # Eliminate any nextSquaresToSearch that don't exist as they're out of the grid.
                realNextSquaresToSearch = [nsq for nsq in nextSquaresToSearch if (nsq[0] >= 0 and nsq[0] < squares.shape[0] and nsq[1] >= 0 and nsq[1] < squares.shape[1])]

                squaresToSearch = realNextSquaresToSearch

                ringIdxRad += 1

            # Now we've found at least one possible tri, find the one with the closest
            # centroid (probably not the best approach but simple and should work
            # fine for now).
            centroidDists = np.array([np.linalg.norm(triangles[t].centroid - point) for t in possibleTris])
            triToUse = possibleTris[np.argmin(centroidDists)]

            return triangles[triToUse].calc_interpolated_value(point)

        assert(False), 'Error: get_piecewise_linear_interpolation_function_for_tri_mesh has gone wrong somehow and a node has slipped through and not been assigned a triangle.'
        return

    # Now return the mapping function we've generated.
    return interpolating_function






#############################################################################
#############################################################################


def get_triangle_finder_function(nodes, triangulation):
    # You feed this function a set of 2D nodes and a triangulation that define a triangle
    # mesh. This function returns a function that can be fed a 2D point [x,y] 
    # and outputs the tri in the mesh that that point lies in, or a sensible tri
    # close to the point if the point lies outside the mesh (though the point should
    # be only slightly outside the mesh for this to work well).
    # (This function is almost a copy-paste of get_piecewise_linear_interpolation_function_for_tri_mesh.)
    # I expect this function could be made a lot faster, possibly by removing some
    # of the class structure. So a standard usage might look like
    # triangle_finder = create_triangle_finder_function(nodes, triangulation)
    # for point in list_of_xy_points_that_are_not_nodes_of_the_mesh:
    #   triangle_id_corresponding_to_this_point = triangle_finder(point)


    num_nodes = nodes.shape[0]
    num_tris = triangulation.shape[0]

    # There may or may not be a slightly faster though less intuitive way
    # to do this using barycentric coords - see various answers on
    # Stackoverflow 2049582, where I also found the method used here.
    # It works as follows: pick any vertex and draw the two sides outwards
    # from it. For the point to be in the triangle, the the corresponding
    # 'cross-products' must have opposite sign. Then look at the sign of
    # the final side's cross product with a vector from a different vertex
    # to the point to finish the determination.
    # I think in rare cases where a point is almost exactly on an edge,
    # this method can fail (see
    # http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html)
    # and you'll have to do more work in that case. For now I just check
    # that there is no such problem i.e. all nodes are assigned to at least
    # one triangle. Something you could do to detect the edge case is catch
    # when an edge think a point in on its left (say) when going along the
    # edge in either direction. Then you could pick some convention etc to
    # handle it.
    def isPointInTriangle(point, vert0, vert1, vert2):

        vert0ToPoint_x = point[0] - vert0[0]
        vert0ToPoint_y = point[1] - vert0[1]

        check1 = (vert1[0]-vert0[0]) * vert0ToPoint_y - (vert1[1]-vert0[1]) * vert0ToPoint_x > 0

        if ( (vert2[0]-vert0[0]) * vert0ToPoint_y - (vert2[1]-vert0[1]) * vert0ToPoint_x > 0 ) == check1:

            return False

        if ( (vert2[0]-vert1[0]) * (point[1] - vert1[1]) - (vert2[1]-vert1[1]) * (point[0] - vert1[0]) > 0 ) != check1:

            return False

        return True

    class Tri:
       def  __init__(self, vertex_labels):
           self.vertex_labels = vertex_labels
           self.vertex_positions = -4321*np.ones((3,2), dtype = 'float')
           self.edges = 4321*np.ones((3,2), dtype = 'float')
           self.centroid = -12345*np.ones(2)
           return

       def set_vertex_positions(self, vertex_positions):
           self.vertex_positions = vertex_positions
           return
       
       # Reorder the storing of vertex labels and positions so the increasing
       # index from 0 to 2 corresponds to traversing the tri anticlockwise.
       def reorder_vertices(self):
           
           vec_0_to_1 = self.vertex_positions[1,:] -self.vertex_positions[0,:]
           vec_0_to_2 = self.vertex_positions[2,:] -self.vertex_positions[0,:]
           
           if vec_0_to_1[0]*vec_0_to_2[1] - vec_0_to_1[1]*vec_0_to_2[0] < 0:
               
               temp = self.vertex_labels[1]
               self.vertex_labels[1] = self.vertex_labels[2]
               self.vertex_labels[2] = temp
               
               temp2 = self.vertex_positions[1,:]
               self.vertex_positions[1,:] = self.vertex_positions[2,:]
               self.vertex_positions[2,:] = temp2
           

       def calc_edges(self):
           # The first edge is the vector pointing from vertex 0 to vertex 1. 
           # The second points from vertex 1 to vertex 2. The third points from
           # vertex 2 to vertex 0. We calculate the edges after putting the 
           # vertices in anticlockwise order, so the calculated edges will all
           # also point anticlockwise. Each ROW of the edges matrix is an edge vector.
           self.edges[0,:] = self.vertex_positions[1,:] - self.vertex_positions[0,:]
           self.edges[1,:] = self.vertex_positions[2,:] - self.vertex_positions[1,:]
           self.edges[2,:] = self.vertex_positions[0,:] - self.vertex_positions[2,:]


       def set_centroid(self, centroid):
           self.centroid = centroid
           return



    triangles = np.empty(num_tris, dtype=object)
    for t in range(0, num_tris):
        triangles[t] = Tri(triangulation[t,:])
        triangles[t].set_vertex_positions(nodes[triangulation[t,:],:])
        triangles[t].reorder_vertices()
        triangles[t].calc_edges()
        triangles[t].set_centroid( np.sum(triangles[t].vertex_positions, axis = 0) / 3)


    min_edge_length = np.linalg.norm(triangles[0].edges[0,:])
    for t in range(0, num_tris):
        length = np.linalg.norm(triangles[t].edges[0,:])
        if min_edge_length > length:
            min_edge_length = length
        length = np.linalg.norm(triangles[t].edges[1,:])
        if min_edge_length > length:
            min_edge_length = length
        length = np.linalg.norm(triangles[t].edges[2,:])
        if min_edge_length > length:
            min_edge_length = length



    gridSquareSize = 1.1 * min_edge_length # This I'm sure is only roughly optimal at best, and could be tweaked.
    # The 2*min_edge_length's are largely cautionary, and may be unnecessary overkill.
    grid_xmin = np.amin(nodes[:,0]) - 2 * min_edge_length
    grid_xmax_approx = np.amax(nodes[:,0]) + 2 * min_edge_length
    gridPoint_xvals = np.linspace( grid_xmin, grid_xmin + ceil((grid_xmax_approx-grid_xmin)/gridSquareSize)*gridSquareSize, ceil((grid_xmax_approx-grid_xmin)/gridSquareSize)+1 )
    grid_xmax = np.amax(gridPoint_xvals)

    grid_ymin = np.amin(nodes[:,1]) - 2 * min_edge_length
    grid_ymax_approx = np.amax(nodes[:,1]) + 2 * min_edge_length
    gridPoint_yvals = np.linspace( grid_ymin, grid_ymin + ceil((grid_ymax_approx-grid_ymin)/gridSquareSize)*gridSquareSize, ceil((grid_ymax_approx-grid_ymin)/gridSquareSize)+1 )
    grid_ymax = np.amax(gridPoint_yvals)


    xx, yy = np.meshgrid(gridPoint_xvals, gridPoint_yvals, sparse=False, indexing='ij') # Matrices of x and y grid point positions respectively



    class Square:
       def  __init__(self, vertPoints):
           self.Idxs = []
           self.vertPoints = vertPoints
           sides = np.zeros((2,4))
           for s in range(0,4):
               sides[:,s] = vertPoints[(s+1)%4,:] - vertPoints[s,:]
           self.sides = sides
           self.triIntersections = []
           return

       def set_Idxs(self, Idxs):
           self.Idxs = Idxs
           return

       def append_TriIntersection(self, triLabel):
           if (not (triLabel in self.triIntersections)): # if statement is precautionary, should never be necessary in this code.
               self.triIntersections.append(triLabel)
           return

    squares = np.empty((xx.shape[0]-1,xx.shape[1]-1), dtype=object)
    for i in range(0, squares.shape[0]):
        for j in range(0, squares.shape[1]):
            squares[i,j] = Square(np.array([[xx[i,j], yy[i,j]], [xx[i+1,j], yy[i+1,j]], [xx[i+1,j+1], yy[i+1,j+1]], [xx[i,j+1], yy[i,j+1]]]))
            squares[i,j].set_Idxs([i,j])
    numSquares = squares.size
    print('Number of squares in grid for get_piecewise_linear_interpolation_function_for_tri_mesh is ' +str(numSquares))

    def SquareIdxsOfPoint(point):
        # Returns the i,j indexes in the squares array corresponding to the square
        # that a given x,y point is in.
        # assert(point[0] >= grid_xmin and point[1] >= grid_ymin and point[0] <= grid_xmax and point[1] <= grid_ymax), 'The point you provide to SquareIdxsOfPoint() must lie within the grid of squares!'
        if( point[0] >= grid_xmin and point[0] <= grid_xmax  and point[1] >= grid_ymin and point[1] <= grid_ymax ): # if point lies within the grid of squares
            return (floor((point[0]-grid_xmin) / gridSquareSize),  floor((point[1]-grid_ymin) / gridSquareSize))
        else:
            return False # point does not lie within the grid of squares


    def convexPolyIntersectionCheck(points1, edges1, points2, edges2):
        # Function that returns bool indicating whether two convex 2D polygons
        # intersect/overlap, based on separating axis theorem. Each ROW of the
        # points1,2 inputs should hold point coords, going round anticlockwise, and
        # the nth COLUMN of edges1,2 should be a side vector from the nth point to
        # the n+1th (so also anticlockwise).
        # See Ericson 'Real-Time Collision Detection'.
        # Some optimisation of this might be possible, see e.g. http://web.archive.org/web/20141127210836/http://content.gpwiki.org/index.php/Polygon_Collision

        # I have not thought about edge cases here where the polygons
        # exactly share a point or nearly do and floating point causes problems etc,
        # as it really shouldn't matter for the current application.
        pointSepVecs = np.empty((points1.shape[0], points2.shape[0], 2))
        for p1 in range(0, points1.shape[0]):
            for p2 in range(0, points2.shape[0]):
                    pointSepVecs[p1,p2,:] = points2[p2,:] - points1[p1,:]


        # Check if any of polygon 1's edges give a separating axis
        for e in range(0, edges1.shape[1]):

            abandonThisEdge = False

            for p in range(0, points2.shape[0]):

                # For an edge to give a separating axis, all the other polygon's
                # points must lie to the right of the edge as you going along it
                # according to our anticlockwise convention.
                # If this edge does not give a separating axis, move on.
                if int(copysign(1, pointSepVecs[e,p,0] * edges1[1,e] - pointSepVecs[e,p,1] * edges1[0,e])) == -1:

                    abandonThisEdge = True
                    break

            if abandonThisEdge == True:
                continue
            # If a separating axis is found, the polygons do not intersect.
            return False

         # Check if any of polygon 2's edges give a separating axis
        for e in range(0, edges2.shape[1]):

            abandonThisEdge = False

            for p in range(0, points1.shape[0]):

                # If this edge does not give a separating axis, move on.
                if int(copysign(1, pointSepVecs[p,e,0] * edges2[1,e] - pointSepVecs[p,e,1] * edges2[0,e])) == 1:

                    abandonThisEdge = True
                    break

            if abandonThisEdge == True:
                continue

            return False

        # If no separating axis found, the polygons intersect
        return True

    # Let the squares know which triangles intersect them.
    for t in range(0, num_tris):
        min_x_vert = triangles[t].vertex_positions[ np.argmin(triangles[t].vertex_positions[:,0]), :]
        max_x_vert = triangles[t].vertex_positions[ np.argmax(triangles[t].vertex_positions[:,0]), :]
        min_y_vert = triangles[t].vertex_positions[ np.argmin(triangles[t].vertex_positions[:,1]), :]
        max_y_vert = triangles[t].vertex_positions[ np.argmax(triangles[t].vertex_positions[:,1]), :]

        # We test all squares in the grid-aligned rectangle defined by the tri verts
        idxsToCheck_x =  range(SquareIdxsOfPoint(min_x_vert)[0], SquareIdxsOfPoint(max_x_vert)[0]+1)
        idxsToCheck_y =  range(SquareIdxsOfPoint(min_y_vert)[1], SquareIdxsOfPoint(max_y_vert)[1]+1)

        for i in idxsToCheck_x:
            for j in idxsToCheck_y:

                if convexPolyIntersectionCheck(triangles[t].vertex_positions, triangles[t].edges.transpose(), squares[i,j].vertPoints, squares[i,j].sides) == True:
                    squares[i,j].append_TriIntersection(t)



    def triangle_finder(point):
        squarePointIsInIfAny = SquareIdxsOfPoint(point)
        # If point is in a square...
        if squarePointIsInIfAny != False:

            # If point's square intersects any triangles...
            possibleTris = squares[squarePointIsInIfAny].triIntersections

            for t in possibleTris:
                # If point is in THIS triangle, return that triangle's ID.
                if isPointInTriangle(point,  triangles[t].vertex_positions[0,:], triangles[t].vertex_positions[1,:], triangles[t].vertex_positions[2,:]) == True:

                    return t

        # Otherwise (i.e. the square had no triangle
        # intersections, or it did but the point wasn't in any of those triangles),
        # we instead just find the closest square to the point that does have
        # intersections, and use the triangle it intersects whose centroid is
        # closest to the point. You could be more careful and really find e.g. the
        # closest point in the set of triangles to the given point, and use the
        # corresponding triangle. However the approach here should be fine for
        # our application. For now we exclude cases where the point fell outside
        # the square grid entirely, but that should be fixed in future, by imagining
        # the grid did extend out to where the point is. The expanding-ring-of-squares
        # search is also done a bit wastefully- in particular when it reaches the edge
        # of the grid in one direction it could just stop rather than checking new
        # 'possibleSquare's in that direction for potential adding to nextSquaresToSearch.
        if squarePointIsInIfAny == False:
                assert(False), 'Error: For now we need all points you want to map to be within the square grid, so make the grid bigger if need be.'
        else:
            possibleTris = []
            squaresToSearch = [squarePointIsInIfAny]
            startSquare = squarePointIsInIfAny
            ringIdxRad = 0 # Search based on a square 'ring' of squares, with the outer squares this many indices from the central start square

            while True:

                for sq in squaresToSearch:
                    possibleTris.extend(squares[sq].triIntersections)

                if len(possibleTris) != 0:
                    break

                nextSquaresToSearch = []
                for sq in squaresToSearch:

                    if sq[0] == startSquare[0] - ringIdxRad:
                        nextSquaresToSearch.append((sq[0]-1, sq[1]))

                    elif sq[0] == startSquare[0] + ringIdxRad:
                        nextSquaresToSearch.append((sq[0]+1, sq[1]))

                    elif sq[1] == startSquare[1] - ringIdxRad:
                        nextSquaresToSearch.append((sq[0], sq[1]-1))

                    elif sq[1] == startSquare[1] + ringIdxRad:
                        nextSquaresToSearch.append((sq[0], sq[1]+1))

                # Corners of ring add an extra square
                nextSquaresToSearch.append((startSquare[0] + ringIdxRad + 1, startSquare[1] + ringIdxRad + 1))
                nextSquaresToSearch.append((startSquare[0] - ringIdxRad - 1, startSquare[1] + ringIdxRad + 1))
                nextSquaresToSearch.append((startSquare[0] - ringIdxRad - 1, startSquare[1] - ringIdxRad - 1))
                nextSquaresToSearch.append((startSquare[0] + ringIdxRad + 1, startSquare[1] - ringIdxRad - 1))

                # Eliminate any nextSquaresToSearch that don't exist as they're out of the grid.
                realNextSquaresToSearch = [nsq for nsq in nextSquaresToSearch if (nsq[0] >= 0 and nsq[0] < squares.shape[0] and nsq[1] >= 0 and nsq[1] < squares.shape[1])]

                squaresToSearch = realNextSquaresToSearch

                ringIdxRad += 1

            # Now we've found at least one possible tri, find the one with the closest
            # centroid (probably not the best approach but simple and should work
            # fine for now).
            centroidDists = np.array([np.linalg.norm(triangles[t].centroid - point) for t in possibleTris])
            triToUse = possibleTris[np.argmin(centroidDists)]

            return triToUse

        assert(False), 'Error: get_triangle_finder_function has gone wrong somehow because a point has slipped through not been assigned a triangle.'
        return

    # Now return the mapping function we've generated.
    return triangle_finder


##############################################################################
##############################################################################

def write_Line_VTK(fileName, headerText, points, dotted=False):
    # Function to take an Nx3 array of point positions representing a set of 
    # straight line segments joined end to end like point0 -- point1 -- point2 --...
    # etc, and write these to an output file, in a (subset of) the .vtk legacy 
    # format for reading into into ParaView. This is useful for drawing a line 
    # in the reference state, mapping it to the final state using e.g.
    # get_map2DPointToSurface_Function, and then drawing that deformed curve 
    # onto the deformed surface in ParaView.
    # If the optional argument dotted == True, then the line segments are 
    # instead taken to be point0 -- point1, point2 -- point3 etc, to create a 
    # dotted line.

    # Set desired format and precision of floating point quantities.
    fp_fmt = '%.7e'

    numPoints = points.shape[0]

    # Some simple checks on a few things. These are not exhaustive.
    assert(
    isinstance(fileName, str)
    and isinstance(headerText, str)
    and points.shape[1] == 3
    ), "Error: Something is wrong with one of the input arguments to write_Line_VTK. Aborting."


    with open(fileName, 'w') as file:
        
        #Write preamble and point coordinates to file.
        preamble = "# vtk DataFile Version 2.0\n" + headerText + "\n" + "ASCII\n" + "DATASET POLYDATA\n" + "POINTS %d double" % numPoints
        np.savetxt(file, points, fmt=fp_fmt, delimiter=' ', header=preamble, comments='')

        if dotted == False:
    
            # Now do line connectivities with  preamble. First column of '2's is 
            # part of format, indicating that each row is a line segment connecting
            # just two points. The overall line is the union of all these segments.
            preamble = "LINES %d %d" % (numPoints-1, 3*(numPoints-1))
            np.savetxt(file, np.column_stack((2*np.ones(numPoints-1, dtype=int), np.arange(0, numPoints-1), np.arange(1, numPoints))), fmt='%d', delimiter=' ', header=preamble, comments='')
            
        elif dotted == True:
            
            assert(numPoints % 2 == 0), "Error: To draw a dotted line you must supply an even number of points."

            preamble = "LINES %d %d" % (numPoints//2, 3*(numPoints//2))
            np.savetxt(file, np.column_stack(( 2*np.ones(numPoints//2, dtype=int), np.arange(0, numPoints, 2), np.arange(0, numPoints, 2)+1 )), fmt='%d', delimiter=' ', header=preamble, comments='')

    return




##############################################################################
##############################################################################



def multiArrayElementwiseExtremum(min_or_max_string, arrayTuple):
    # Function that takes a string that's either `min' or 'max', and a tuple of
    # numpy arrays, and perform np.minimum or np.maximum on the first two input 
    # arrays,  and and then on the result of that and the third input array ...
    # .... etc, so you end up with the element-wise min or max of the whole 
    # tuple of arrays.
    
    assert( 
    min_or_max_string == 'min' or  min_or_max_string == 'max'
    ), "Error: Something is wrong with one of the input arguments to multiArrayElementwiseExtremum. Aborting."
    
    numArrays = len(arrayTuple)
    
    output = arrayTuple[0]
    
    for a in range(1, numArrays):
        
        if( min_or_max_string == 'min' ):
            output = np.minimum(output, arrayTuple[a])
        elif( min_or_max_string == 'max' ):
            output = np.maximum(output, arrayTuple[a])    

    return output


##############################################################################
##############################################################################



def multi_array_elementwise_and_or_or(and_or_or_string, array_tuple):
    # Function that takes a string that's either 'and' or 'or', and a tuple of
    # numpy arrays, and perform np.logical_and or np.logical_or on the first two input 
    # arrays,  and and then on the result of that and the third input array ...
    # .... etc, so you end up with the element-wise min or max of the whole 
    # tuple of arrays.
    
    assert( 
    and_or_or_string == 'and' or  and_or_or_string == 'or'
    ), "Error: Something is wrong with one of the input arguments to multi_array_elementwise_and_or_or. Aborting."
    
    num_arrays = len(array_tuple)
    
    output = array_tuple[0]
    
    for a in range(1, num_arrays):
        
        if( and_or_or_string == 'and' ):
            output = np.logical_and(output, array_tuple[a])
        elif( and_or_or_string == 'or' ):
            output = np.logical_or(output, array_tuple[a])    

    return output
