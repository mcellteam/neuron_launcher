import bpy, bmesh
from mathutils import Vector
import math

import sys

# Main

def f_close_open_caps(context):

    print("> Running: f_close_open_caps")

    # Get the active object
    ob_list = context.selected_objects
    
    if len(ob_list) == 1:
        ob = ob_list[0]

        # List of verts, edges, faces to add
        v_add_list = []
        e_add_list = []
        f_add_list = []

        # Select all edges
        bpy.ops.object.mode_set(mode='EDIT')
        context.tool_settings.mesh_select_mode = (False, True, False)
        bpy.ops.mesh.select_all(action='SELECT')

        # Select boundary loop
        bpy.ops.mesh.region_to_loop()
        # Nobody knows why but need to flip here
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')

        # All boundaries
        all_bdry_edges = []
        for i_e, ed in enumerate(ob.data.edges):
            if ed.select == True:
                all_bdry_edges.append(i_e)

        '''
        print("---")
        print("All boundary edges")
        print(all_bdry_edges)
        '''

        # Go through each boundary and select one
        while len(all_bdry_edges) > 0:
            
            print("Checking a loop - length remaining: " + str(len(all_bdry_edges)))

            # Make sure to start by deselect all
            bpy.ops.mesh.select_all(action='DESELECT')

            # Select one
            bpy.ops.object.mode_set(mode='OBJECT')
            ob.data.edges[all_bdry_edges[0]].select = True
            bpy.ops.object.mode_set(mode='EDIT')
        
            # Get connected edge loop
            bpy.ops.mesh.loop_multi_select(ring=False)
            # For whatever reason, need to flip object/edit after this
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.mode_set(mode='EDIT')

            # What is selected
            edge_loop = []
            for i_e, edge in enumerate(ob.data.edges):
                if edge.select == True:
                    edge_loop.append(i_e)

            # Center
            vert_loop = [list(ob.data.edges[item].vertices) for item in edge_loop]
            vert_u_loop = list(set([item for sublist in vert_loop for item in sublist]))
            vert_co_loop = [ob.data.vertices[item].co for item in vert_u_loop]
            ctr = vert_co_loop[0]
            for v_co in vert_co_loop[1:]:
                ctr += v_co
            ctr /= len(vert_co_loop)

            # Convert center to world space (right now it is object)
            # smwi = ob.matrix_world.inverted()
            ctr = ob.matrix_world * ctr # World space

            # Construct connect vertex list
            vert_conn_loop = [vert_loop[0][0], vert_loop[0][1]]
            del vert_loop[0]
            while len(vert_loop) > 0:
                # Find connecting piece
                for i_vs, vs in enumerate(vert_loop):
                    if vs[0] == vert_conn_loop[-1]:
                        vert_conn_loop.append(vs[1])
                        del vert_loop[i_vs]
                        break
                    if vs[1] == vert_conn_loop[-1]:
                        vert_conn_loop.append(vs[0])
                        del vert_loop[i_vs]
                        break

            # Add center to the vert list
            v_add_list.append(ctr)
            ctr_id = len(ob.data.vertices) + len(v_add_list) - 1

            # Add edges, faces
            # First face/edge
            f_add_list.append([ctr_id, vert_conn_loop[0], vert_conn_loop[1]])
            e_add_list.append([ctr_id, vert_conn_loop[0]])
            e_add_list.append([ctr_id, vert_conn_loop[1]])
            # Other faces/edges
            for i in range(1,len(vert_conn_loop)-1):
                f_add_list.append([ctr_id, vert_conn_loop[i],vert_conn_loop[i+1]])
                e_add_list.append([ctr_id, vert_conn_loop[i+1]])
                               
            # Delete from all_bdry_edges
            edge_loop.sort()
            edge_loop.reverse()
            for e in edge_loop:
                del all_bdry_edges[all_bdry_edges.index(e)]

            print("Finished this loop - length remaining: " + str(len(all_bdry_edges)))

        print("Finished all loops - updating object")

        # Make the vertex, edge, face lists
        v_list = [item.co for item in ob.data.vertices] + v_add_list
        v_list = [tuple(item) for item in v_list]
        e_list = [item.vertices for item in ob.data.edges] + e_add_list
        f_list = [item.vertices for item in ob.data.polygons] + f_add_list

        # Make a new obj
        mesh_new = bpy.data.meshes.new(ob.name+"_closed_mesh")
        mesh_new.from_pydata(v_list,e_list,f_list)
        mesh_new.validate(verbose=False) # Important! and i dont know why
        mesh_new.update()
        obj_new = bpy.data.objects.new(ob.name+"_closed",mesh_new)
        context.scene.objects.link(obj_new)

    print("> Finished: f_close_open_caps")






