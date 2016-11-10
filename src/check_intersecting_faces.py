import bpy, bmesh
from mathutils import Vector
import math

import sys

def intersect3D_RayTriangle(r, t):
    u = t[1] - t[0]
    v = t[2] - t[0]
    n = u.cross(v)
    if (n == Vector([0,0,0])):
        return -1 # Degenerate

    direc = r[1] - r[0]
    w0 = r[0] - t[0]
    a = -1.0*n.dot(w0)
    b = n.dot(direc)
    if (abs(b) < 0.00000000001): # ray is parallel to triangle
        if (a==0): # lies in triangle plane
            return 2
        else: # degenerate
            return 0

    rl = a / b
    if rl < 0.0 or rl > 1.0:
        return 0

    ipt = r[0] + rl * direc

    uu = u.dot(u)
    uv = u.dot(v)
    vv = v.dot(v)
    w = ipt - t[0]
    wu = w.dot(u)
    wv = w.dot(v)
    dist = uv * uv - uu * vv

    s = (uv * wv - vv * wu) / dist
    if (s < 0.0 or s > 1.0):
        return 0 # ipt is outside of triangle
    t = (uv * wu - uu * wv) / dist
    if (t < 0.0 or (s+t) > 1.0):
        return 0 # degenerate

    return 1 # intersect

def check_triangles(t1, t2):
    # Construct all the rays in t1
    for i,j in [(0,1),(0,2),(1,2)]:
        val = intersect3D_RayTriangle([t1[i],t1[j]],t2)
        if val != 0:
            print(val)

# Main

if __name__ == "__main__":

    print("> Running")
    
    # Get the objects
    ob_surf = bpy.data.objects['P40_surface']
    ob_seg = bpy.data.objects['P40_segment']

    # Get the centers of all faces
    v_surf_triplets = [tuple(item.vertices) for item in ob_surf.data.polygons]
    v_surf_co = [[ob_surf.data.vertices[item[0]].co,ob_surf.data.vertices[item[1]].co,ob_surf.data.vertices[item[2]].co] for item in v_surf_triplets]
    v_seg_triplets = [tuple(item.vertices) for item in ob_seg.data.polygons]
    v_seg_co = [[ob_seg.data.vertices[item[0]].co,ob_seg.data.vertices[item[1]].co,ob_seg.data.vertices[item[2]].co] for item in v_seg_triplets]
   
    # Check if any face in surf list matches seg list
    l = len(v_surf_co)
    for i_1, tri_1 in enumerate(v_surf_co):
        if i_1 % 100 == 0:
            print(str(i_1) + " / " + str(l))

        # Check all other centers
        for i_2, tri_2 in enumerate(v_seg_co):
            check_triangles(tri_1, tri_2)

    print("> Finished")


