import bpy, bmesh
from mathutils import Vector
import math

import sys

# Main

if __name__ == "__main__":

    print("> Running")
    
    # Get the objects
    ob_surf = bpy.data.objects['P40_surface']
    ob_seg = bpy.data.objects['P40_segment']

    # Get the centers of all faces
    ctr_surf_list = [item.center for item in ob_surf.data.polygons]
    ctr_seg_list = [item.center for item in ob_seg.data.polygons]

    # Check if any face in surf list matches seg list
    l = len(ctr_surf_list)
    for i,ctr in enumerate(ctr_surf_list):
        if i % 1000 == 0:
            print(str(i) + " / " + str(l))

        # Check all other centers
        for i_o,ctr_o in enumerate(ctr_seg_list):
            len = (ctr-ctr_o).length
            if len < 0.0000001:
                print("These centers match - do the triangles?")
                print(ctr)
                print(ctr_o)
                vs = ob_surf.data.polygons[i]
                for j in range(0,3):
                    print(ob_surf.data.vertices[j].co)
                vs = ob_seg.data.polygons[i_o]
                for j in range(0,3):
                    print(ob_seg.data.vertices[j].co)

    print("> Finished")


