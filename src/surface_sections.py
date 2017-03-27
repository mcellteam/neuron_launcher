import bpy, bmesh
from mathutils import Vector
import math

import collections

import numpy as np

import functools

import itertools

import os

import sys

# To add objects to MCell
from cellblender.cellblender_utils import preserve_selection_use_operator

# Function to project point onto line
'''
def project_pt_line(v, w, p):
    v_to_p = p-v
    v_to_w = w-v
    
    vtw2 = v_to_w.length
    vtp2 = v_to_p.length
    
    v_to_p.normalize()
    v_to_w.normalize()
    
    vtp_dot_vtw = v_to_p.dot(v_to_w)
    
    t2 = vtp_dot_vtw * vtp2 / vtw2;
    
    v_to_w = w-v
    
    return v + (v_to_w * t2)
'''

# Deep copy a list
def unshared_copy(inList):
    if isinstance(inList, list):
        return list( map(unshared_copy, inList) )
    return inList

# Read in the swc file to know what's connected
def get_connections(fname):
    
    # Store the swc data
    swc_data = []
    
    # Read the file
    f = open(fname,'r')
    for i_line,line in enumerate(f):
        line_s = line.split();
        if len(line_s) > 0 and line_s[0] != "#":
            line_data = []
            line_data.append(int(float(line_s[0])))
            line_data.append(int(float(line_s[1])))
            line_data.append(float(line_s[2]))
            line_data.append(float(line_s[3]))
            line_data.append(float(line_s[4]))
            line_data.append(float(line_s[5]))
            line_data.append(int(float(line_s[6])))
            swc_data.append(line_data)

    # Vertex coordinates
    vert_list = [Vector(item[2:5]) for item in swc_data]

    # Radius list
    r_list = [item[5] for item in swc_data]

    # Find connections
    pt_connect = []
    for i in range(0,len(swc_data)):
        pt_connect.append([])

    for i,data in enumerate(swc_data):
        if i > 0: # The first is the -1 index
            pt_connect[data[6]-1].append(i+1)

    pt_connect_dup = unshared_copy(pt_connect) # copy

    for i_list, sub_list in enumerate(pt_connect_dup):
        for idx in sub_list:
            pt_connect[idx-1].append(i_list+1)

    return pt_connect, vert_list, r_list

# Function to construct a neighbour list from a list of face vertices
'''
def construct_neighbor_list(face_list):
    neigh_list = []
    for f in face_list:
        neigh_list.append([])

    # Go through each face
    for i_f, f in enumerate(face_list):
        
        # Go through every edge pair of vertices
        for i in range(0,3):
            # If we already have three neighbors, stop
            if len(neigh_list[i_f]) == 3:
                break
            
            # Edge vertices
            ev1 = f[i]
            ev2 = f[(i+1)%3]
            # Search tediously for this face's neighbor
            for i_fo, fo in enumerate(face_list[i_f+1:]): # and don't double count
                if ev1 in fo and ev2 in fo:
                    # Found the neighbor!
                    neigh_list[i_f].append(i_fo)
                    neigh_list[i_fo].append(i_f)
                    break

    return neigh_list
'''

# Construct planes for vertices
def construct_dividing_plane_normals(pt_connect, vert_co_list):

    # Dictionary of vertex to segment to normals
    # Each index is a vertex in the SWC
    # Each entry is another dictionary: one for each segment that connects to this vertex
    # Each entry in this dictionary is a list of all the plane normals to consider
    normal_dict = {}
    
    # Go through all vertices
    for i_pt, conn_pts in enumerate(pt_connect):
        pt = i_pt + 1
        normal_dict[pt] = {}
    
        # Go through all connecting points
        for i_pt_o, pt_o in enumerate(conn_pts[:-1]):

            # Vector for this segment
            vec_1 = vert_co_list[pt_o-1] - vert_co_list[pt-1]
            vec_1.normalize()
            
            # Check against all other connecting segments
            for pt_p in conn_pts[i_pt_o+1:]:

                # Vector for this segment
                vec_2 = vert_co_list[pt_p-1] - vert_co_list[pt-1]
                vec_2.normalize()

                # Average
                vec_a = 0.5*(vec_1 + vec_2)

                # Cross
                vec_c = vec_1.cross(vec_2)

                # Plane's normal
                vec_n = vec_a.cross(vec_c)

                # Ensure that we are pointing in the same hemisphere as vec_1 or 2
                if vec_n.dot(vec_1) < 0:
                    vec_n *= -1.0

                # Store
                if pt_o in normal_dict[pt]:
                    normal_dict[pt][pt_o].append(vec_n)
                else:
                    normal_dict[pt][pt_o] = [vec_n]
                if pt_p in normal_dict[pt]:
                    normal_dict[pt][pt_p].append(-1.0*vec_n)
                else:
                    normal_dict[pt][pt_p] = [-1.0*vec_n]

                # Make a plane for viz
                '''
                vec_a.normalize()
                vec_c.normalize()
                r = 1.0
                verts = [vert_co_list[pt-1]+r*vec_a, vert_co_list[pt-1]+r*vec_c, vert_co_list[pt-1]-r*vec_a, vert_co_list[pt-1]-r*vec_c];
                # Make a new obj
                mesh_new = bpy.data.meshes.new("planes_%02d_%02d_%02d_mesh"%(pt,pt_o,pt_p))
                verts = [tuple(item) for item in verts]
                edges = [[0,1],[1,2],[2,3],[3,0]]
                faces = [[0,1,2,3]]
                mesh_new.from_pydata(verts,edges,faces)
                mesh_new.validate(verbose=False) # Important! and i dont know why
                mesh_new.update()
                obj_new = bpy.data.objects.new("planes_%02d_%02d_%02d"%(pt,pt_o,pt_p),mesh_new)
                bpy.context.scene.objects.link(obj_new)
                '''

    return normal_dict

# Reconstruct the connectivity allowed
# Returns a dictionary with keys of sections in the form of tuples (v1 v2) 
# and values of a list of sections that are allowed to border it
def construct_connectivity(pt_connect):

    conn_dict = {}

    for pt1_0, conn_pts in enumerate(pt_connect):
        pt1 = pt1_0 + 1
        for pt2 in conn_pts:
            key = (min(pt1,pt2),max(pt1,pt2))
            conn_dict[key] = []

            # Get the allowed conns
            for i_pt in [pt1,pt2]:
                for j_pt in pt_connect[i_pt-1]:
                    conn_sc = (min(i_pt,j_pt),max(i_pt,j_pt))
                    if not conn_sc in conn_dict[key]:
                        conn_dict[key].append(conn_sc)

    return conn_dict

# Main

def f_surface_sections(context, swc_filepath):

    print("> Running: f_surface_sections")

    # Get the object
    ob_list = context.selected_objects

    if len(ob_list) != 1:
        raise TypeError("Please select only one object.")

    ob = ob_list[0]

    # Get data from the SWC file
    pt_connect, vert_co_list, vert_r_list = get_connections(swc_filepath)
    
    # Construct the connectivity allowed
    conn_dict = construct_connectivity(pt_connect)

    # Construct dividing planes
    normal_dict = construct_dividing_plane_normals(pt_connect, vert_co_list)

    # List of vertices that are endpoints
    endpoint_list = []
    for i_pt_0, conn_pts in enumerate(pt_connect):
        i_pt = i_pt_0 + 1
        if len(conn_pts) == 1:
            endpoint_list.append(i_pt)

    # List of sections
    sc_list = []
    # Go through all the points
    for i_pt, conn_pts in enumerate(pt_connect):
        for conn_pt in conn_pts:
            # Check against duplicates
            if conn_pt > i_pt + 1:
                sc_list.append((i_pt+1,conn_pt))

    # ORDER the sections by average radius of the two vertices, in decreasing order
    # That is, assign faces to large radius sections first, as these deserve "more" faces
    sc_r_list = []
    for v_pair in sc_list:
        r_ave = 0.5*(vert_r_list[v_pair[0]-1] + vert_r_list[v_pair[1]-1])
        sc_r_list.append(r_ave)

    # Zipped sort
    sc_r_list, sc_list = (list(t) for t in zip(*sorted(zip(sc_r_list, sc_list),reverse=True)))

    # Create a list of face indexes in ob.data.polygons that need to be assigned to sections
    # As faces are assigned they are removed from this list
    face_idx_list = list(range(0,len(ob.data.polygons)))

    # zero vector
    zv = Vector([0.0,0.0,0.0])

    # Store the face idxs assigned to each section
    sc_face_dict = {}

    # Store the vertices that make up the border loop around each section
    sc_brdr_vert_vict = {}

    # Go through every section and assign faces to it
    for i_ctr,sc in enumerate(sc_list):
        print("Assigning faces for section " + str(sc) + " (" + str(i_ctr+1) + "/" + str(len(sc_list)) + ")...")
        
        # Face selection mode
        context.tool_settings.mesh_select_mode = (False, False, True)

        ### 
        # Step 1: Disregard faces that violate our bordering condition from the SWC file
        ###

        # Disregard list of faces that include these vertices
        vert_disregard_list = []

        # Go through all existing sections
        for sc_existing, vert_brdr_list in sc_brdr_vert_vict.items():

            # Is our current section allowed to border this one?
            if not sc_existing in conn_dict[sc]:

                # No, it isn't allowed to border this one!

                # Store to disregard
                vert_disregard_list += vert_brdr_list

        # Remove duplicates (is this necessary or superfluous?)
        vert_disregard_list = list(set(vert_disregard_list))

        ###
        # Step 2: find faces that are possible candidates for belonging to this section
        ###

        # Coordinate of vertices of this section
        v_co_0 = vert_co_list[sc[0]-1]
        v_co_1 = vert_co_list[sc[1]-1]

        # Max dist that faces are allowed to be away from either vertex
        max_dist = 1.5*(v_co_1-v_co_0).length

        # Get all the planes from this section, that occur at this vertex
        plane_normals_st = []
        for i_vert in [0,1]:
            normal_tmp = normal_dict[sc[i_vert]]
            if not sc[(i_vert+1)%2] in normal_tmp:
                # This vertex has no other sections that connect to it!
                # Skip and check the next section
                plane_normals_st.append([])
            else:
                plane_normals_st.append(normal_tmp[sc[(i_vert+1)%2]])

        # First find all faces that are within max_dist
        tri_list = []
        for f in face_idx_list:
    
            # Distances from each vert
            ctr = ob.data.polygons[f].center

            disps = []
            if not sc[1] in endpoint_list:
                disps.append(ctr - v_co_0)
            else:
                disps.append(zv) # Free pass for endpoint
            if not sc[0] in endpoint_list:
                disps.append(ctr - v_co_1)
            else:
                disps.append(zv) # Free pass for endpoint

            # Check that distances are allowed
            if disps[0].length <= max_dist and disps[1].length <= max_dist:
                    
                # Does this face share vertices with any of the illegal vertices? If so, disregard
                vert_f_idxs = list(ob.data.polygons[f].vertices)
                if not vert_f_idxs[0] in vert_disregard_list and not vert_f_idxs[1] in vert_disregard_list and not vert_f_idxs[2] in vert_disregard_list:

                    # This is not sufficient - this just indicates that this face lies in the intersection region of two spheres!
                    # Now use the normals of the planes that define the regions to further constrain which elements are allowed
                    
                    # Flag
                    COUNT_FACE = True
                    
                    # Check each vertex
                    for i_vert in [0,1]:

                        # But don't check the endpoints
                        if not sc[i_vert] in endpoint_list:
                        
                            # Check all the normals
                            for normal in plane_normals_st[i_vert]:
                        
                                if normal.dot(disps[i_vert]) < 0:
                                    COUNT_FACE = False
                                    break
                    
                            # Don't keep search for this face if we already violated one of the vertices
                            # i.e. break out of the i_vert = [0,1] loop
                            if COUNT_FACE == False:
                                break
                                    
                    # Append
                    if COUNT_FACE == True:
                        tri_list.append(f)

        ###
        # Step 3: Proceed using the "select more" = bpy.ops.mesh.select_more() function
        ###

        # First set the mode
        bpy.ops.object.mode_set(mode='EDIT')
        # Face selection mode
        context.tool_settings.mesh_select_mode = (False, False, True)
        # Deselect all
        bpy.ops.mesh.select_all(action='DESELECT')
        # Must be in object mode to select faces!
        bpy.ops.object.mode_set(mode='OBJECT')

        # Initially select the faces in tri_list which are closest to each of the vertices
        min_f0 = -1
        dist_f0 = -1
        min_f1 = -1
        dist_f1 = -1
        for f in tri_list:
            dist_v0 = (ob.data.polygons[f].center - v_co_0).length
            if dist_v0 < dist_f0 or dist_f0 == -1:
                min_f0 = f
                dist_f0 = dist_v0
            dist_v1 = (ob.data.polygons[f].center - v_co_1).length
            if dist_v1 < dist_f1 or dist_f1 == -1:
                min_f1 = f
                dist_f1 = dist_v1

        ob.data.polygons[min_f0].select = True
        ob.data.polygons[min_f1].select = True

        # Search for faces that belong to this section

        # Init list of faces that belong to this section
        sc_face_dict[sc] = []
        # Init list of indexes to delete from the face_idx_list to prevent double checking
        delete_list = []

        SEARCH = True
        while SEARCH:
            
            # Use the select more function
            # Must be in edit mode to use select more!
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_more()
            bpy.ops.object.mode_set(mode='OBJECT')
            
            # Check each of the possible faces if it is selected
            new_faces = []
            delete_tri_list = []
            
            # Fast but worse
            '''
            for i_f,f in enumerate(tri_list):
                if ob.data.polygons[f].select == True:
            '''
            # Slower but better
            for f, f_face in enumerate(ob.data.polygons):
                if f_face.select == True:
                    if f in tri_list:
                        i_f = tri_list.index(f)
                        # Add it as a valid face
                        new_faces.append(f)
                        # Store for deletion from tri_list to prevent double counting at next iteration
                        delete_tri_list.append(i_f)
                    else:
                        # Turn off the face selection so that we do not continue to "grow" in this direction using the select more function
                        f_face.select = False

            # Check that there was something new added
            if len(new_faces) == 0:
                SEARCH = False
                break
            
            # Delete triangles from tri_list to prevent double counting in the future (+ faster!)
            delete_tri_list.reverse()
            for i_f in delete_tri_list:
                del tri_list[i_f]

            # Add the new faces
            for f in new_faces:
                sc_face_dict[sc].append(f)
                delete_list.append(face_idx_list.index(f))

        # Select the region
        '''
        # First set the mode
        bpy.ops.object.mode_set(mode='EDIT')
        # Face selection mode
        context.tool_settings.mesh_select_mode = (False, False, True)
        # Deselect all
        bpy.ops.mesh.select_all(action='DESELECT')
        # Must be in object mode to select faces!
        bpy.ops.object.mode_set(mode='OBJECT')
        # Select
        for t in sc_face_dict[sc]:
            ob.data.polygons[t].select = True
        bpy.ops.object.mode_set(mode='EDIT')
        '''

        # Do the deletion
        delete_list.sort()
        delete_list.reverse()
        for i_f in delete_list:
            del face_idx_list[i_f]

        ###
        # Step 4 - Get the vertices that make up the border of this section
        ###

        # First deselect all
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        # Next select the faces that belong to this section
        for f in ob.data.polygons:
            if f.index in sc_face_dict[sc]:
                f.select = True

        # Convert to edge loop
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.region_to_loop()

        # Get the edges
        bpy.ops.object.mode_set(mode='OBJECT')
        edge_loop = [i.index for i in ob.data.edges if i.select]

        # Convert to vertices and store
        sc_brdr_vert_vict[sc] = [item for edge in edge_loop for item in ob.data.edges[edge].vertices]


    ###
    # Finished assigning faces to sections!
    ###

    ###
    # Sanity check: are all faces assigned?
    ###

    f_ctr = 0
    for f_list in sc_face_dict.values():
        f_ctr += len(f_list)

    print("Number of faces to assign: " + str(len(ob.data.polygons)))
    print("Number of faces assigned: " + str(f_ctr))

    if len(ob.data.polygons) != f_ctr:
        print("WARNING! There are faces that have not been assigned to regions.")

    ###
    # Turn the dictionary into MCell regions
    ###

    print("Adding regions to MCell...")

    # Ensure object is active
    context.scene.objects.active = ob
    ob.select = True
    
    # First add the object to MCell
    preserve_selection_use_operator(bpy.ops.mcell.model_objects_add, ob)

    # Add the sections as regions
    existing_reg_names = [item.name for item in ob.mcell.regions.region_list]
    for sc_id,f_list in sc_face_dict.items():

        # Region name
        reg_name = "sc_%02d_%02d"%sc_id

        # Check that the region does not exist
        if reg_name in existing_reg_names:
            # Get region
            new_reg = ob.mcell.regions.region_list[reg_name]
            # Clear it
            new_reg.reset_region(context)
        else:
            # Make region
            ob.mcell.regions.add_region_by_name(context,reg_name)
            # Get region
            new_reg = ob.mcell.regions.region_list[reg_name]
    
        # Assign faces
        new_reg.set_region_faces(ob.data, f_list)

    # Update (but dont validate because who knows)
    ob.data.update()

    print("> Finished: f_surface_sections")


def f_surface_sections_V2(context, swc_filepath):

    print("> Running: f_surface_sections")

    # Get the object
    ob_list = context.selected_objects

    if len(ob_list) != 1:
        raise TypeError("Please select only one object.")

    ob = ob_list[0]

    # Get data from the SWC file
    pt_connect, vert_co_list, vert_r_list = get_connections(swc_filepath)
    
    # Construct the connectivity allowed
    conn_dict = construct_connectivity(pt_connect)

    # Construct dividing planes
    normal_dict = construct_dividing_plane_normals(pt_connect, vert_co_list)

    # List of vertices that are endpoints
    endpoint_list = []
    for i_pt_0, conn_pts in enumerate(pt_connect):
        i_pt = i_pt_0 + 1
        if len(conn_pts) == 1:
            endpoint_list.append(i_pt)

    # List of sections
    sc_list = []
    # Go through all the points
    for i_pt, conn_pts in enumerate(pt_connect):
        for conn_pt in conn_pts:
            # Check against duplicates
            if conn_pt > i_pt + 1:
                sc_list.append((i_pt+1,conn_pt))

    # ORDER the sections by average radius of the two vertices, in decreasing order
    # That is, assign faces to large radius sections first, as these deserve "more" faces
    sc_r_list = []
    for v_pair in sc_list:
        r_ave = 0.5*(vert_r_list[v_pair[0]-1] + vert_r_list[v_pair[1]-1])
        sc_r_list.append(r_ave)

    # Zipped sort
    sc_r_list, sc_list = (list(t) for t in zip(*sorted(zip(sc_r_list, sc_list),reverse=True)))

    # Store the face idxs assigned to each section
    sc_face_dict = {}

    # Store the vertices that make up the border loop around each section
    sc_brdr_vert_vict = {}

    # Go through every section and assign faces to it
    for i_ctr,sc in enumerate(sc_list):
        print("Assigning faces for section " + str(sc) + " (" + str(i_ctr+1) + "/" + str(len(sc_list)) + ")...")
        
        # Face selection mode
        context.tool_settings.mesh_select_mode = (False, False, True)

        ### 
        # Step 1: Disregard faces that violate our bordering condition from the SWC file
        ###

        # Disregard list of faces that include these vertices
        vert_disregard_list = []

        # Go through all existing sections
        for sc_existing, vert_brdr_list in sc_brdr_vert_vict.items():

            # Is our current section allowed to border this one?
            if not sc_existing in conn_dict[sc]:

                # No, it isn't allowed to border this one!

                # Store to disregard
                vert_disregard_list += vert_brdr_list

        # Remove duplicates (is this necessary or superfluous?)
        vert_disregard_list = list(set(vert_disregard_list))

        ###
        # Step 2: find faces that are possible candidates for belonging to this section
        ###

        # Coordinate of vertices of this section
        v_co_0 = vert_co_list[sc[0]-1]
        v_co_1 = vert_co_list[sc[1]-1]

        # Max dist that faces are allowed to be away from either vertex
        max_dist = 1.5*(v_co_1-v_co_0).length

        # Get all the planes from this section, that occur at this vertex
        plane_normals_st = []
        for i_vert in [0,1]:
            normal_tmp = normal_dict[sc[i_vert]]
            if not sc[(i_vert+1)%2] in normal_tmp:
                # This vertex has no other sections that connect to it!
                # Skip and check the next section
                plane_normals_st.append([])
            else:
                plane_normals_st.append(normal_tmp[sc[(i_vert+1)%2]])

        # First find all faces that are within max_dist
        tri_list = []
        for f in face_idx_list:
    
            # Distances from each vert
            ctr = ob.data.polygons[f].center

            disps = []
            if not sc[1] in endpoint_list:
                disps.append(ctr - v_co_0)
            else:
                disps.append(zv) # Free pass for endpoint
            if not sc[0] in endpoint_list:
                disps.append(ctr - v_co_1)
            else:
                disps.append(zv) # Free pass for endpoint

            # Check that distances are allowed
            if disps[0].length <= max_dist and disps[1].length <= max_dist:
                    
                # Does this face share vertices with any of the illegal vertices? If so, disregard
                vert_f_idxs = list(ob.data.polygons[f].vertices)
                if not vert_f_idxs[0] in vert_disregard_list and not vert_f_idxs[1] in vert_disregard_list and not vert_f_idxs[2] in vert_disregard_list:

                    # This is not sufficient - this just indicates that this face lies in the intersection region of two spheres!
                    # Now use the normals of the planes that define the regions to further constrain which elements are allowed
                    
                    # Flag
                    COUNT_FACE = True
                    
                    # Check each vertex
                    for i_vert in [0,1]:

                        # But don't check the endpoints
                        if not sc[i_vert] in endpoint_list:
                        
                            # Check all the normals
                            for normal in plane_normals_st[i_vert]:
                        
                                if normal.dot(disps[i_vert]) < 0:
                                    COUNT_FACE = False
                                    break
                    
                            # Don't keep search for this face if we already violated one of the vertices
                            # i.e. break out of the i_vert = [0,1] loop
                            if COUNT_FACE == False:
                                break
                                    
                    # Append
                    if COUNT_FACE == True:
                        tri_list.append(f)

        ###
        # Step 3: Proceed using the "select more" = bpy.ops.mesh.select_more() function
        ###

        # First set the mode
        bpy.ops.object.mode_set(mode='EDIT')
        # Face selection mode
        context.tool_settings.mesh_select_mode = (False, False, True)
        # Deselect all
        bpy.ops.mesh.select_all(action='DESELECT')
        # Must be in object mode to select faces!
        bpy.ops.object.mode_set(mode='OBJECT')

        # Initially select the faces in tri_list which are closest to each of the vertices
        min_f0 = -1
        dist_f0 = -1
        min_f1 = -1
        dist_f1 = -1
        for f in tri_list:
            dist_v0 = (ob.data.polygons[f].center - v_co_0).length
            if dist_v0 < dist_f0 or dist_f0 == -1:
                min_f0 = f
                dist_f0 = dist_v0
            dist_v1 = (ob.data.polygons[f].center - v_co_1).length
            if dist_v1 < dist_f1 or dist_f1 == -1:
                min_f1 = f
                dist_f1 = dist_v1

        ob.data.polygons[min_f0].select = True
        ob.data.polygons[min_f1].select = True

        # Search for faces that belong to this section

        # Init list of faces that belong to this section
        sc_face_dict[sc] = []
        # Init list of indexes to delete from the face_idx_list to prevent double checking
        delete_list = []

        SEARCH = True
        while SEARCH:
            
            # Use the select more function
            # Must be in edit mode to use select more!
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_more()
            bpy.ops.object.mode_set(mode='OBJECT')
            
            # Check each of the possible faces if it is selected
            new_faces = []
            delete_tri_list = []
            
            # Fast but worse
            '''
            for i_f,f in enumerate(tri_list):
                if ob.data.polygons[f].select == True:
            '''
            # Slower but better
            for f, f_face in enumerate(ob.data.polygons):
                if f_face.select == True:
                    if f in tri_list:
                        i_f = tri_list.index(f)
                        # Add it as a valid face
                        new_faces.append(f)
                        # Store for deletion from tri_list to prevent double counting at next iteration
                        delete_tri_list.append(i_f)
                    else:
                        # Turn off the face selection so that we do not continue to "grow" in this direction using the select more function
                        f_face.select = False

            # Check that there was something new added
            if len(new_faces) == 0:
                SEARCH = False
                break
            
            # Delete triangles from tri_list to prevent double counting in the future (+ faster!)
            delete_tri_list.reverse()
            for i_f in delete_tri_list:
                del tri_list[i_f]

            # Add the new faces
            for f in new_faces:
                sc_face_dict[sc].append(f)
                delete_list.append(face_idx_list.index(f))

        # Select the region
        '''
        # First set the mode
        bpy.ops.object.mode_set(mode='EDIT')
        # Face selection mode
        context.tool_settings.mesh_select_mode = (False, False, True)
        # Deselect all
        bpy.ops.mesh.select_all(action='DESELECT')
        # Must be in object mode to select faces!
        bpy.ops.object.mode_set(mode='OBJECT')
        # Select
        for t in sc_face_dict[sc]:
            ob.data.polygons[t].select = True
        bpy.ops.object.mode_set(mode='EDIT')
        '''

        # Do the deletion
        delete_list.sort()
        delete_list.reverse()
        for i_f in delete_list:
            del face_idx_list[i_f]

        ###
        # Step 4 - Get the vertices that make up the border of this section
        ###

        # First deselect all
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.object.mode_set(mode='OBJECT')

        # Next select the faces that belong to this section
        for f in ob.data.polygons:
            if f.index in sc_face_dict[sc]:
                f.select = True

        # Convert to edge loop
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.region_to_loop()

        # Get the edges
        bpy.ops.object.mode_set(mode='OBJECT')
        edge_loop = [i.index for i in ob.data.edges if i.select]

        # Convert to vertices and store
        sc_brdr_vert_vict[sc] = [item for edge in edge_loop for item in ob.data.edges[edge].vertices]


    ###
    # Finished assigning faces to sections!
    ###

    ###
    # Sanity check: are all faces assigned?
    ###

    f_ctr = 0
    for f_list in sc_face_dict.values():
        f_ctr += len(f_list)

    print("Number of faces to assign: " + str(len(ob.data.polygons)))
    print("Number of faces assigned: " + str(f_ctr))

    if len(ob.data.polygons) != f_ctr:
        print("WARNING! There are faces that have not been assigned to regions.")

    ###
    # Turn the dictionary into MCell regions
    ###

    print("Adding regions to MCell...")

    # Ensure object is active
    context.scene.objects.active = ob
    ob.select = True
    
    # First add the object to MCell
    preserve_selection_use_operator(bpy.ops.mcell.model_objects_add, ob)

    # Add the sections as regions
    existing_reg_names = [item.name for item in ob.mcell.regions.region_list]
    for sc_id,f_list in sc_face_dict.items():

        # Region name
        reg_name = "sc_%02d_%02d"%sc_id

        # Check that the region does not exist
        if reg_name in existing_reg_names:
            # Get region
            new_reg = ob.mcell.regions.region_list[reg_name]
            # Clear it
            new_reg.reset_region(context)
        else:
            # Make region
            ob.mcell.regions.add_region_by_name(context,reg_name)
            # Get region
            new_reg = ob.mcell.regions.region_list[reg_name]
    
        # Assign faces
        new_reg.set_region_faces(ob.data, f_list)

    # Update (but dont validate because who knows)
    ob.data.update()

    print("> Finished: f_surface_sections")


