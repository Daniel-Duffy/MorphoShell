#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: dld34


"""

import numpy as np
import distmesh as dm
from math import sqrt, atan2, pi, sin, tan, cos, asin, cosh, sinh, log, asinh, floor, acos, atan
import matplotlib.pyplot as plt
from scipy import optimize

from mesh_processing_module import (
    write_VTK,
    prog_Tensors_From_Top_And_Bottom_Director_Patterns,
    read_Between_Lines,
    get_map2DPointToSurface_Function,
    multiArrayElementwiseExtremum,
    optimize_mesh_tri_and_node_labels)




# Set working directory. Can be relative to python script.
#working_dir = '/home/dlduffy/c++_stuff/'


# Radius of initial flat disk.
R = 1.0

   
# Define signed distance function used to make mesh.
dist_func = lambda p: dm.dcircle(p, 0, 0, R)

# Uniform is probably not bad for now.
coarsening = lambda p: dm.huniform(p)

# Sets the shortest edge length in the mesh
approx_shortest_edge_length = R / 50

fixed_points = None

# Create mesh using PyDistMesh. Nodes holds node positions. The triangulation 
# gives node labels (row indices in 'nodes') starting from 0.
nodes, triangulation = dm.distmesh2d(dist_func, coarsening, approx_shortest_edge_length, (-1.1*R, -1.1*R, 1.1*R, 1.1*R), fixed_points)

print('Finished creating mesh.')
optimize_mesh_tri_and_node_labels(triangulation, nodes) # Makes the tri and node labels correspond to their spatial locations to some extent, which can improve performance.
num_nodes = nodes.shape[0]
num_tris = triangulation.shape[0]
print('Number of nodes = '+str(num_nodes))
print('Number of triangles = '+str(num_tris))

# Get array of triangle centroid coordinates.
tri_centroids = np.zeros((num_tris,2))
for i in range(0, num_tris):
    tri_centroids[i,:] = ( nodes[triangulation[i,0],:] + nodes[triangulation[i,1],:] + nodes[triangulation[i,2],:] ) / 3
    


# Setting these to negative numbers means the values are set to a uniform single number chosen in the settings file instead of here.
ref_thicknesses = -1234.0 * np.ones(num_tris)
ref_shear_moduli = -1234.0 * np.ones(num_tris)

# Constraint indicators and other tags.
constraint_indicators = np.zeros(num_nodes)
node_tags = np.zeros(num_nodes)

tri_tags = np.zeros(num_tris)

# Standard Cartesian and polar coordinates defined in the flat 2D plane of the initial sheet.
x = tri_centroids[:,0]
y = tri_centroids[:,1]
r = np.sqrt(x**2 + y**2)
phi = np.arctan2(y, x) + pi # Between 0 and 2pi.


# Angle between in-plane director the x-axis, so in the flat state director = [cos(director_angle), sin(director_angle)].
director_angles = phi + np.pi/2 # Azimuthal director.

# Lambda is the contraction factor along the in-plane director. Lambda is dialled linearly
# in time from 1 down to Lambda_final.
Lambda_final = 0.85

# Contraction factor in direction perpendicular to director is Lambda**(-nu_therm) 
nu_therm = 0.0

# Create some arrays; we'll fill them next.
abar_info = np.zeros((num_tris,3))
bbar_info = np.zeros((num_tris,3))

abar_info = np.stack((director_angles, Lambda_final*np.ones(num_tris), nu_therm*np.ones(num_tris)), axis=-1)
# We'll send the above info about abar to MorphoShell (c++), and within that c++ simulation code the following
# situation will be arranged (but we don't need to actually do it in python, which is why I've commented it out).
# It's the standard protocol when the contraction direction lies in the plane of the material. 
"""
for t in range(0, num_tris):
	abar_diag = np.array([[Lambda**2, 0],
                             [0, Lambda**(-2*nu_therm)]])
	rot_mat = np.array([[np.cos(director_angles[t]), -np.sin(director_angles[t])], 
                           [np.sin(director_angles[t]),  np.cos(director_angles[t])]])
	abar = rot_mat @ abar_diag @ (rot_mat.T)
	abar_comps[t,:] = np.array([abar[0,0],abar[0,1],abar[1,1]])
"""

# Now we do bbar. We'll take the preferred principal curvature directions to lie parallel to the director and 
# perpendicular to the director. This might (for example) correspond to a setup in which the top of the sheet contracts by more 
# than the bottom of the sheet (but along the same director on both sides). Thus there are two preferred principal curvatures, kappabar_par 
# and kappabar_perp. In the c++ code, if using the dialling protocol labelled "4" (and with dialling_from_ansatz_rather_than_ref_state = false 
# in the settings file), both kappabar_par and kappabar_perp will  be dialled linearly over time, from 0 to their final values defined below.
for t in range(0, num_tris):    

	kappabar_par_final = 2.5 / R

	kappabar_perp_final = 0.0 / R

 	bbar_info[t,:] = np.array([director_angles[t], kappabar_par_final, kappabar_perp_final])


        
write_VTK('../files_from_python/cone.vtk', 
    ' ', 
    nodes, 
    triangulation, 
    abar_info, 
    bbar_info, 
    ref_thicknesses, 
    ref_shear_moduli,
    tri_tags,
    constraint_indicators, 
    node_tags)


# Create ansatz file.
ansatz_nodes = np.zeros((num_nodes,3))
r = np.sqrt(nodes[:,0]**2 + nodes[:,1]**2)
for i in range(0, num_nodes):
    
    # Slightly curved state.
    ansatz_nodes[i,0] = nodes[i,0]
    ansatz_nodes[i,1] = nodes[i,1]
    ansatz_nodes[i,2] = -(0.001/R) * (nodes[i,0]**2 + nodes[i,1]**2)
    
    
   
write_VTK('../files_from_python/cone_ansatz.vtk', 
    'dial_factor = 1.0', 
    ansatz_nodes, 
    triangulation, 
    abar_info=np.zeros((num_tris,3)), 
    bbar_info=np.zeros((num_tris,3)), 
    ref_thicknesses=np.zeros(num_tris), 
    ref_shear_moduli=np.zeros(num_tris), 
    tri_tags=np.zeros(num_tris), 
    constraint_indicators=np.zeros(num_nodes), 
    node_tags=np.zeros(num_nodes)  )








