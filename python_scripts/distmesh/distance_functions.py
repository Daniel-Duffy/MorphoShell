# encoding: utf-8
"""Distance functions."""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np

__all__ = [
    # Distance functions:
    'dsector',
    'dblock',
    'dcircle',
    'ddiff',
    #'dellipse',
    #'dellipsoid',
#    'dexpr',
    'dintersect',
    'dmatrix3d',
    'dmatrix',
    'dpoly', # Duffy uncommented this
    'drectangle',
    'drectangle0',
    'dsegment', # Duffy uncommented this
    'dsphere',
    'dunion',

    # Mesh size functions:
    'hmatrix3d',
    'hmatrix',
    'huniform',

    # Generic node manipulation:
    'protate',
    'pshift',
    
    # Other
    'clamp',
    ]

# These are used very often:
min = np.minimum
max = np.maximum

def clamp(x, minval, maxval):
    return np.minimum( np.maximum(x, minval), maxval)

#-----------------------------------------------------------------------------
# Signed distance functions
#-----------------------------------------------------------------------------

# Distance to a sector of a circle, or what you might call a pizza slice.
# Added by Daniel Duffy, essentially copying the amazing Inigo Quilez https://iquilezles.org/articles/distfunctions2d/
# The argument ang is the half-angle of the sector, while rad is the radius.
def dsector(p, ang, rad):
    q = np.stack((np.fabs(p[:,0]),p[:,1]), axis=-1)
    c = np.array([np.sin(ang), np.cos(ang)])
    l = np.sqrt((q**2).sum(-1)) - rad
    temp = clamp(q@c, 0.0, rad)
    temp2 = np.stack((temp*c[0], temp*c[1]), axis=-1)
    m = np.sqrt((( q - temp2 )**2).sum(-1))
    return np.maximum(l, m*np.sign(q[:,0]*c[1]-q[:,1]*c[0]))
    
def dblock(p,x1,x2,y1,y2,z1,z2):
    return -min(min(min(min(min(-z1+p[:,2],z2-p[:,2]),-y1+p[:,1]),y2-p[:,1]),-x1+p[:,0]),x2-p[:,0])

def dcircle(p,xc,yc,r):
    """Signed distance to circle centered at xc, yc with radius r."""
    return np.sqrt(((p-np.array([xc,yc]))**2).sum(-1))-r

def ddiff(d1,d2):
    """Signed distance to set difference between two regions described by
    signed distance functions d1 and d2.

    Not exact the true signed distance function for the difference,
    for example around corners.
    """
    return max(d1,-d2)

#from distmesh._distance_functions import dellipse

#from distmesh._distance_functions import dellipsoid

# dexpr not implemented

def dintersect(d1,d2):
    """Signed distance to set intersection of two regions described by signed
    distance functions d1 and d2.

    Not exact the true signed distance function for the intersection,
    for example around corners.
    """
    return max(d1,d2)

def dmatrix3d(p,xx,yy,zz,dd):
    """Signed distance function by interpolation of the values dd on the
    Cartesian grid xx, yy, zz."""
    return ml.interp3_linear(xx,yy,zz,dd,p[:,0],p[:,1],p[:,2])

def dmatrix(p,xx,yy,dd):
    """Signed distance function by interpolation of the values dd on the
    Cartesian grid xx, yy."""
    return ml.interp2_linear(xx,yy,dd,p[:,0],p[:,1])

def dpoly(p,pv):
    """Signed distance function for polygon with vertices pv.

    Usually pv should also be provided as fixed points in distmesh2d.

    pv should be provided as a list of coordinates [(x0,y0), (x1,y1), ...]
    or an array of shape (nv, 2).
    
    Added by Duffy: When defining a closed polygon, you should repeat a vertex,
    so the last row of pv should be the same as the first.
    """
    from matplotlib.path import Path
    return (-1)**Path(pv).contains_points(p) * dsegment(p, pv).min(1)

def drectangle0(p,x1,x2,y1,y2):
    """Signed distance function for rectangle with corners (x1,y1), (x2,y1),
    (x1,y2), (x2,y2).

    See drectangle for a simpler version ignoring corners.
    """
    d1=y1-p[:,1]
    d2=-y2+p[:,1]
    d3=x1-p[:,0]
    d4=-x2+p[:,0]

    d5=np.sqrt(d1**2+d3**2)
    d6=np.sqrt(d1**2+d4**2)
    d7=np.sqrt(d2**2+d3**2)
    d8=np.sqrt(d2**2+d4**2)

    d=-min(min(min(-d1,-d2),-d3),-d4)

    ix=(d1>0)*(d3>0)
    d[ix]=d5[ix]
    ix=(d1>0)*(d4>0)
    d[ix]=d6[ix]
    ix=(d2>0)*(d3>0)
    d[ix]=d7[ix]
    ix=(d2>0)*(d4>0)
    d[ix]=d8[ix]

    return d

def drectangle(p,x1,x2,y1,y2):
    """Signed distance function for rectangle with corners (x1,y1), (x2,y1),
    (x1,y2), (x2,y2).

    This has an incorrect distance to the four corners. See drectangle0 for a
    true distance function.
    """
    return -min(min(min(-y1+p[:,1],y2-p[:,1]),-x1+p[:,0]),x2-p[:,0])

#from distmesh._distance_functions import dsegment
# Copied and translated by Daniel Duffy from /home/daniel/helpful_stuff/temp/pydistmesh-master/distmesh/src/distance_functions.c
# Input format: p is a num_points x 2 array as usual, while pv is a num_vertices x 2 array of vertex positions.
# The first line segment has endpoints pv[0,:] and pv[1,:], while the second has endpoints pv[1,:] and pv[2,:], etc.
def dsegment(p, pv):
    
    dists_to_segments = np.zeros((p.shape[0], pv.shape[0]-1))
    
    # Loop over segments. There's probably a more pythonic and faster way than 
    # using a loop, using broadcasting or maybe np.add.outer or maybe even einsum. 
    # See e.g. 
    # https://stackoverflow.com/a/33848814
    # https://stackoverflow.com/questions/70370522
    # I doubt any such approach is really worth doing for the numbers of segments
    # we're likely to use though.
    for s in range(0, pv.shape[0]-1):
        
        p1x = pv[s,0]
        p1y = pv[s,1]
        p2x = pv[s+1,0]
        p2y = pv[s+1,1]
    
        c1 = (p2x-p1x)*(p[:,0]-p1x) + (p2y-p1y)*(p[:,1]-p1y) # dot2(v,w)
        c2 = (p2x-p1x)*(p2x-p1x) + (p2y-p1y)*(p2y-p1y) # dot2(v,v)
            
        c1_greaterthan_0 = c1 > 0
        selector = np.logical_not(c1_greaterthan_0)
        dists_to_segments[selector,s] = np.sqrt((p[selector,0]-p1x)**2 + (p[selector,1]-p1y)**2)
        selector = np.logical_and(c1_greaterthan_0, c1>=c2)
        dists_to_segments[selector,s] = np.sqrt((p[selector,0]-p2x)**2 + (p[selector,1]-p2y)**2)
        selector = np.logical_and(c1_greaterthan_0, c1<c2)
        dists_to_segments[selector,s] = np.sqrt((p[selector,0]-(p1x+c1[selector]/c2*(p2x-p1x)))**2 + (p[selector,1]-(p1y+c1[selector]/c2*(p2y-p1y)))**2)
            
    return dists_to_segments


def dsphere(p,xc,yc,zc,r):
    """Signed distance function for a sphere centered at xc,yc,zc with radius
    r."""
    return np.sqrt((p[:,0]-xc)**2+(p[:,1]-yc)**2+(p[:,2]-zc)**2)-r

def dunion(d1,d2):
    """Signed stance function for the set union of two regions described by
    signed distance functions d1, d2.

    This not a true signed distance function for the union, for example around
    corners.
    """
    return min(d1,d2)

#-----------------------------------------------------------------------------
# Mesh size functions
#-----------------------------------------------------------------------------

def hmatrix3d(p,xx,yy,zz,dd,hh):
    """Mesh size function by interpolation of the values hh on the Cartesian
    grid xx, yy, zz."""
    return ml.interp3_linear(xx,yy,zz,hh,p[:,0],p[:,1],p[:,2]);

def hmatrix(p,xx,yy,dd,hh):
    """Mesh size function by interpolation of the values hh on the Cartesian
    grid xx, yy."""
    return ml.interp2_linear(xx,yy,hh,p[:,1],p[:,2]);

def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])

#-----------------------------------------------------------------------------
# Generic node manipulation
#-----------------------------------------------------------------------------

def protate(p,phi):
    """Rotate points p the angle phi around origin."""
    A = np.array(((np.cos(phi), -np.sin(phi)),
                  (np.sin(phi),  np.cos(phi))))
    return p.dot(A)

def pshift(p,x0,y0):
    """Move points p by (x0,y0)."""
    return p + [x0,y0]
