"""


"""

import numpy as np

from math import sqrt, atan2, acos, sin, cos, pi

import matplotlib.pyplot as plt

import distmesh as dm

from mesh_processing_module import (
    write_VTK,
    prog_Tensors_From_Top_And_Bottom_Director_Patterns,
    read_Between_Lines,
    get_map2DPointToSurface_Function)

from scipy import optimize
    

#######################################



# Set working directory. Can be relative to python script.
workingDir = '/home/daniel/c++_stuff/gauss-codazzi-mainardi/files_from_python'



Ri = 1.0
Ro = 1.4

   
# Define signed distance function used to make mesh.
dist_func = lambda p: dm.ddiff( dm.dcircle(p, 0, 0, Ro), dm.dcircle(p, 0, 0, Ri) )

# Annulus with an break in it
dist_func = lambda p: dm.ddiff(
                            dm.ddiff( dm.dcircle(p, 0, 0, Ro), dm.dcircle(p, 0, 0, Ri) ),
                            dm.dpoly(p, np.array([[0.99*Ri, 0],
                                                  [1.01*Ro, 0],
                                                  [1.01*Ro, 0.01*Ro],
                                                  [0.99*Ri, 0.01*Ro],
                                                  [0.99*Ri, 0]]))
                            )


coarsening = lambda p: dm.huniform(p)
#coarsening = lambda p:  dm.dcircle(p, 0, 0, 0) / Rcyl

approx_shortest_edge_length = (Ro-Ri) / 20
   

# Create mesh using PyDistMesh. Nodes holds node positions. The triangulation
# gives node labels (row indices in 'nodes') starting from 0.
nodes, triangulation = dm.distmesh2d(dist_func, coarsening, approx_shortest_edge_length, (-1.1*Ro,-1.1*Ro,1.1*Ro,1.1*Ro)) 
print('Finished creating mesh.')
num_nodes = nodes.shape[0]
num_tris = triangulation.shape[0]
print('Number of nodes = '+str(num_nodes))
print('Number of triangles = '+str(num_tris))
   
# Get array of triangle centroid coordinates.
tri_centroids = np.zeros((num_tris,2))
for t in range(0, num_tris):
    tri_centroids[t,:] = ( nodes[triangulation[t,0],:] + nodes[triangulation[t,1],:] + nodes[triangulation[t,2],:] ) / 3
  
   
  
    

r = np.sqrt(tri_centroids[:,0]**2 + tri_centroids[:,1]**2)
phi = np.arctan2(tri_centroids[:,1], tri_centroids[:,0]) + pi # Between 0 and 2pi.

abar_info = np.ones((num_tris,3))
abar_info[:,1] = 0

bbar_info = np.zeros((num_tris,3))
for t in range(0, num_tris):

    phihat = np.array([-sin(phi[t]), cos(phi[t])])
    
    princ_curv = 0.5 / r[t]
    
    bbar = princ_curv * np.outer(phihat, phihat)
    
    bbar_info[t,:] = np.array([bbar[0,0],bbar[0,1],bbar[1,1]])

        


ref_thicknesses = -1234.0 * np.ones(num_tris)
ref_shear_moduli = -1234.0 * np.ones(num_tris)
   	
   	
# Constraint indicators and other tags.
constraint_indicators = np.zeros(num_nodes)
node_tags = np.zeros(num_nodes)

tri_tags = np.zeros(num_tris)

   
write_VTK(workingDir+'/annulus_ref.vtk', 
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
ansatz_nodes[:,0] = nodes[:,0] 
ansatz_nodes[:,1] = nodes[:,1] 
ansatz_nodes[:,2] = - 0.02 * r**2 / Ro**2
 
 
write_VTK(workingDir+'/annulus_ansatz.vtk', 
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


print("REMEMBER: ANSATZ SMOOTHING and cuthill mckee IN PYTHON, GENERALIZED ALPHA DYNAMICS IN C++")
print(" store both storage orders of big sparse matrices to aid parallelisation ")
print(" inexact jacobian could be a good idea, Grinspun et al -- Simple and efficient implementation of discrete plates and shells, p.26 first column, and also Grinspun et al -- Cubic Shells.")
"""
People sometimes say that finite element stiffness matrices are 'usually banded', with no justification
for this claim. I'm confident it boils down only to this: suppose your nodes are arranged spatially like
 0  1  2  3  4  5  6  7  8  9
10 11 12 13 14 15 16 17 18 19
20 21 22 23 24 25 26 27 28 29
where the numbers are the node labels/IDs. You can see that if your domain is about n nodes wide 
(left-to-right, n=10 in this example) and m nodes tall, then spatially neighbouring nodes will have labels
differing by around n or less. But the square stiffness matrix is (n*m) x (n*m). So e.g if m=n, the length of 
a row of the stiffness matrix is n^2, and the non-zero elements on that row will only span about n either side 
of the diagonal, so for large n you actually do have a rather banded matrix, even though you numbered your nodes
only in a very naive way!
I'll do scipy's reverse cuthill mckee on the tri-tri adjacency matrix (which I computed nicely in c++, just copy the algo over), to
rename the tris. Then invent a whole new node numbering, so that tri 0's first vertex is node 0, its second is node 1, and so on. You just 
keep a running dictionary of old->new node ids, and for every new node you look at you see whether it already has a new node 
number assigned to its original node number (and switch it to its new one if so), or else you just add a new entry to the dictionary,
with the new node number that is the lowest not yet used.
"""

#######################################################
#######################################################
# Ansatzing using other runs' output files.
ansatzWithOtherRun = False
if( ansatzWithOtherRun == True ):
    
    ref_filename_for_ansatz = '/home/daniel/c++_stuff/mingchao_annuli/making_initial_state/input_files_used/annulus_ref.vtk'
    deformed_filename_for_ansatz = '/home/daniel/c++_stuff/mingchao_annuli/making_initial_state/stepcount_560479_output.vtk'
    
    
    refNodes = read_Between_Lines(ref_filename_for_ansatz, "POINTS", "POLYGONS", return_np_array_using_dtype = 'float')
    refNodes = refNodes[:,0:2]
    
    simTriangulation = read_Between_Lines(ref_filename_for_ansatz, "POLYGONS", "CELL_DATA", return_np_array_using_dtype = 'int')
    simTriangulation = simTriangulation[:,1:]
    
    finalNodes = read_Between_Lines(deformed_filename_for_ansatz, "POINTS", "POLYGONS", return_np_array_using_dtype = 'float')
    
    map2DPointToSurface = get_map2DPointToSurface_Function(refNodes, finalNodes, simTriangulation)
    
    ansatz_nodes = np.zeros((num_nodes,3))
    for n in range(0, num_nodes):
        ansatz_nodes[n,:] = map2DPointToSurface(nodes[n,:])
    
    
    write_VTK(workingDir+'/annulus_ansatz.vtk', 
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
    
    






