import bpy, bmesh
from mathutils import Vector
import math

import collections

import numpy as np

import functools

import itertools

import os

# To add objects to MCell
from cellblender.cellblender_utils import preserve_selection_use_operator

# Time
import time

# Bisect
import bisect

import sys

# Class for a section
class MN_section:
    
    # Init
    def __init__(self, *args, **kwargs):
        # Section name
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = ""
        
        # Tuple of two vertex ids that define this section
        if 'sc_id' in kwargs:
            self.sc_id = kwargs['sc_id']
        else:
            self.sc_id = (-1,-1)
        
        # Boundary cable model indices, points
        if 'sc_pts' in kwargs:
            self.sc_pts = kwargs['sc_pts']
        else:
            self.sc_pts = (Vector([0,0,0]),Vector([0,0,0]))
        
        # Neighboring section tuples at each end pt
        if 'nghbr_sc_ids' in kwargs:
            self.nghbr_sc_ids = kwargs['nghbr_sc_ids']
        else:
            self.nghbr_sc_ids = [[],[]]

        # Number of segments
        self.n_seg = 0
        
        # List of boundary planes between sections at each vertex
        self.planes_sc_brdrs = [[],[]]

        # List of segment boundary planes
        self.planes_sg_brdrs = []

        # List face indexes of segments on the surface, one list for each segment
        self.surf_sg_f_list = []

        # Edge list of vertices
        if 'sc_edge_list' in kwargs:
            self.sc_edge_list = kwargs['sc_edge_list']
        else:
            self.sc_edge_list = []

# Class for a plane
class MN_plane:
    
    # Compute the surface area
    def surf_area(self):
        sa = 0.0
        for f in self.face_list:
            sa += 0.5*( (Vector(self.vert_list[f[1]]) - Vector(self.vert_list[f[0]])).cross(Vector(self.vert_list[f[2]]) - Vector(self.vert_list[f[0]])) ).length
        
        return sa

    '''

    # Make into a separate object
    def make_separate_object(self, obj_vert_list):

        # List of unique vertex ids
        v_all_list = [item for sublist in self.edge_list for item in sublist]
        v_unique_list = list(set(v_all_list))

        # Construct the vertex list
        v_list = []
        v_id_dict = {}
        for v in v_unique_list:
            v_list.append(tuple(obj_vert_list[v]))
            v_id_dict[v] = len(v_list)-1

        v_list.append(tuple(self.ctr_pt))
        ctr_id = len(v_list) - 1

        # Construct the face list
        f_list = []
        for e in self.edge_list:
            f_list.append([v_id_dict[e[0]],v_id_dict[e[1]],ctr_id])

        # Construct the edge list
        e_list = []
        for e in self.edge_list:
            e_list.append([v_id_dict[e[0]],v_id_dict[e[1]]])
        for v in v_unique_list:
            e_list.append([ctr_id,v_id_dict[v]])

        # Make the object!
        mesh_new = bpy.data.meshes.new(self.name + "_mesh")
        mesh_new.from_pydata(v_list,e_list,f_list)
        # Validate and update
        mesh_new.validate(verbose=False) # Important! and i dont know why
        mesh_new.update()
        # Overwrite existing object
        obj_old = bpy.data.objects.get(self.name)
        if obj_old:
            bpy.context.scene.objects.unlink(obj_old)
            bpy.data.objects.remove(obj_old)
        # New one
        obj_new = bpy.data.objects.new(self.name,mesh_new)
        # Link
        bpy.context.scene.objects.link(obj_new)

    '''

    # Init
    def __init__(self, *args, **kwargs):
        # Name
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = ""
        
        # Triplets (two sc + 1 sg) ids for what's on either side, lower first
        if 'sides_ids' in kwargs:
            self.sides_ids = kwargs['sides_ids']
        else:
            self.sides_ids = ((),())

        # Edge list around the loop
        if 'edge_list' in kwargs:
            self.edge_list = kwargs['edge_list']
        else:
            self.edge_list = []

        # Center point of the face
        if 'ctr_pt' in kwargs:
            self.ctr_pt = kwargs['ctr_pt']
        else:
            self.ctr_pt = Vector([0,0,0])

# Function to project point onto line
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

# Deep copy a list
def unshared_copy(inList):
    if isinstance(inList, list):
        return list( map(unshared_copy, inList) )
    return inList

# Read in the swc file to know what's connected
def get_connections(fname):
    global mn_section_dict
    
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

    # Make a new section for all
    for i, sublist in enumerate(pt_connect):
        pt1 = i+1
        for pt2 in sublist:
            if pt2 > pt1:
                # New section
                sec_name = "sc_%02d_%02d"%(pt1,pt2)
                sec_ids = (pt1,pt2)
                sec_pts = (Vector(swc_data[pt1-1][2:5]),Vector(swc_data[pt2-1][2:5]))
                
                nghbrs_min = []
                for conn_pt in pt_connect[pt1-1]:
                    if conn_pt != pt1:
                        nghbrs_min.append((min(pt1,conn_pt),max(pt1,conn_pt)))
                nghbrs_max = []
                for conn_pt in pt_connect[pt2-1]:
                    if conn_pt != pt2:
                        nghbrs_max.append((min(pt2,conn_pt),max(pt2,conn_pt)))

                # Make it
                sec = MN_section(name=sec_name, sc_id=sec_ids, sc_pts=sec_pts, nghbr_sc_ids=(nghbrs_min,nghbrs_max))
                mn_section_dict[sec_ids] = sec

    return


# Main

def f_compartmentize_cyl(context, swc_filepath, n_seg_plen):

    print("> Running: f_compartmentize_cyl")

    # Dictionary of sections
    global mn_section_dict
    mn_section_dict = {}

    # Time
    t_st = []
    t_st.append(time.time())

    # Get data from the SWC file
    get_connections(swc_filepath)

    # Get the active object
    ob_list = context.selected_objects

    if len(ob_list) != 1:
        raise SystemError("Please select one (and only one) object")
    else:
        ob = ob_list[0]

    # Make sure everything is de-selected before we start
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    # Get the regions of the cube object
    reg_list = ob.mcell.regions.region_list

    ###
    # Get all region's boundaries
    ###

    print("> Retrieving border edges...")

    # Go over all sections, get the border edges
    for sec in reg_list:
        # Get the name
        sec_name = sec.name
        
        # Check that it is a section
        if len(sec_name) == 8 and sec_name[0:3] == 'sc_':
            
            # Get the vert indeces
            pt1 = int(sec_name[3:5])
            pt2 = int(sec_name[6:8])
            pt_min = min(pt1,pt2)
            pt_max = max(pt1,pt2)
            
            # Select the faces
            bpy.ops.object.mode_set(mode='EDIT')
            sec.select_region_faces(context)

            # Get the edges
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.region_to_loop()
            bpy.ops.object.mode_set(mode='OBJECT')
            mn_section_dict[(pt_min,pt_max)].sc_edge_list = [i.index for i in ob.data.edges if i.select]

            # Deselect all
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='DESELECT')
            bpy.ops.object.mode_set(mode='OBJECT')

    # Make sure the main guy is not selected
    ob.select = False
    context.scene.objects.active = ob

    # Time
    t_st.append(time.time())
    print("> Retrieved border edges: time: " + str(t_st[-1]-t_st[-2]))

    ###
    # Make SECTION boundaries
    ###

    print("> Making section boundaries...")

    # Go through all of the sections
    for this_sc_id,this_sc in mn_section_dict.items():
        
        # Border
        this_sc_edge = this_sc.sc_edge_list
        
        # Go through both of end verts
        for i_end,nghbr_sc_ids in enumerate(this_sc.nghbr_sc_ids):
            
            # Go through all of the neighboring sections
            for nghbr_sc_id in nghbr_sc_ids:
                
                # Don't double count
                if nghbr_sc_id > this_sc_id:
                    
                    # Border
                    nghbr_sc_edge = mn_section_dict[nghbr_sc_id].sc_edge_list

                    # Get the elements they share
                    edge_share = list(set(this_sc_edge).intersection(set(nghbr_sc_edge)))

                    if len(edge_share) > 0: # Do they share any?
                        
                        # Make a new plane

                        # Get the center at this endpoint
                        ctr_pt = this_sc.sc_pts[i_end]

                        # Get the actual neighboring section
                        nghbr_sc = mn_section_dict[nghbr_sc_id]

                        # Convert the edge id's to actual vertex pairs
                        edge_pair_list = [ob.data.edges[item].vertices for item in edge_share]
                        edge_list = [[int(ids[0]),int(ids[1])] for ids in edge_pair_list]

                        # Make the plane
                        name_plane = this_sc.name + "_B_" + nghbr_sc.name
                        plane = MN_plane(name=name_plane,sides_names=(this_sc.name,nghbr_sc.name), \
                            sides_ids=(this_sc_id,nghbr_sc_id),edge_list=edge_list,ctr_pt=ctr_pt)

                        # Add the plane to both section boundaries
                        this_sc.planes_sc_brdrs[i_end].append(plane)
                        if nghbr_sc_id[0] == this_sc_id[i_end]:
                            nghbr_sc.planes_sc_brdrs[0].append(plane)
                        elif nghbr_sc_id[1] == this_sc_id[i_end]:
                            nghbr_sc.planes_sc_brdrs[1].append(plane)
                        else:
                            print("Something went wrong!")

    # Time
    t_st.append(time.time())
    print("> Made section boundaries: time: " + str(t_st[-1]-t_st[-2]))

    ###
    # Make SEGMENT boundaries
    ###

    print("> Making segment boundaries...")

    # Edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')

    # Go through every section
    for sec in reg_list:
        # Get the name
        sec_name = sec.name
        
        print("Checking section: " + str(sec_name))

        # Check that it is a section
        if len(sec_name) == 8 and sec_name[0:3] == 'sc_':
            
            # Get the sc id
            pt1 = int(sec_name[3:5])
            pt2 = int(sec_name[6:8])
            pt_min = min(pt1,pt2)
            pt_max = max(pt1,pt2)
            this_sc_id = (pt_min,pt_max)

            # Get the sc
            this_sc = mn_section_dict[this_sc_id]

            # Determine the number of segments to make
            sc_unit_vec = this_sc.sc_pts[1] - this_sc.sc_pts[0]
            sc_unit_vec.normalize()
            sc_length = (this_sc.sc_pts[1] - this_sc.sc_pts[0]).length
            n_seg = int(sc_length * n_seg_plen) + 1

            print("Number of segments: " + str(n_seg))

            # Select the faces in this segment
            sec.select_region_faces(context)

            if n_seg != 1:
                    
                # List of dividing lengths from sc_pts[0]
                sc_dividing_len_list = [0.0]
                for i in range(1,n_seg):
                    sc_dividing_len_list.append(i * sc_length / n_seg)
                sc_dividing_len_list.append(2.0*sc_length) # factor 2 = some buffer

                # Make a list of all the faces in this region's projected distances and their compartments
                proj_face_list = []
                for f in ob.data.polygons:
                    if f.select == True:
                        proj_pt = project_pt_line(this_sc.sc_pts[0],this_sc.sc_pts[1],Vector(f.center))
                        proj_dist = (proj_pt - this_sc.sc_pts[0]).length
                        i_sg = bisect.bisect(sc_dividing_len_list, proj_dist) - 1

                        proj_face_list.append((f.index,proj_dist,i_sg))

                # Assign each face to a segment
                
                # List of face id's in each segment
                sg_f_list = []
                for i in range(0,n_seg):
                    sg_f_list.append([])

                # Go through all faces
                for f_index, proj_dist, i_sg in proj_face_list:
                    sg_f_list[i_sg].append(f_index)

                # Prune empty face lists....
                del_list = []
                for i_sg, f_list in enumerate(sg_f_list):
                    if len(f_list) == 0:
                        del_list.append(i_sg)
                del_list.reverse()
                for i in del_list:
                    del sg_f_list[i]

                # Store for this section
                this_sc.surf_sg_f_list = sg_f_list

                # The revised number of segments
                n_seg_rev = len(sg_f_list)

                if n_seg_rev != n_seg:
                    print("Revised number of segments to: " + str(n_seg_rev) + " after assigning faces")

                # Check that there is more than one segment
                if n_seg_rev > 1:

                    # Go through each segment, get the boundary edge
                    sg_bdry_e_list = []
                    for i in range(0,n_seg_rev):
                        sg_bdry_e_list.append([])
                    for i_sg, f_list in enumerate(sg_f_list):

                        # Deselect all faces
                        bpy.ops.object.mode_set(mode='EDIT') # Must be in edit mode
                        bpy.ops.mesh.select_all(action='DESELECT')

                        # Need to be in object mode to select faces
                        bpy.ops.object.mode_set(mode='OBJECT')

                        # Select all faces in this segment
                        for f in f_list:
                            ob.data.polygons[f].select = True

                        # Get the border loop                    
                        bpy.ops.object.mode_set(mode='EDIT') # Must be in edit mode
                        bpy.ops.mesh.region_to_loop()

                        # Need to be in object mode to select edges
                        bpy.ops.object.mode_set(mode='OBJECT')

                        for e in ob.data.edges:
                            if e.select == True:
                                sg_bdry_e_list[i_sg].append(e.index)

                    ###
                    # Determine boundaries and make bounding planes
                    ###
                    for i_sg in range(0,n_seg_rev-1):
                        bdry_list_1 = sg_bdry_e_list[i_sg]
                        bdry_list_2 = sg_bdry_e_list[i_sg+1]

                        # Check for overlap
                        bdry_overlap = list(set(bdry_list_1).intersection(set(bdry_list_2)))

                        if len(bdry_overlap) > 0:

                            # Convert to edge pair list
                            bdry_e_list_tmp = [ob.data.edges[item].vertices for item in bdry_overlap]
                            bdry_e_list = [[int(ids[0]),int(ids[1])] for ids in bdry_e_list_tmp]

                            ###
                            # METHOD 1 for placing center points
                            ###
                            '''
                            # Just place them evenly spaced along the cable model
                            plane_ctr_pt = this_sc.sc_pts[0] + (i_sg+1) * (sc_length/n_seg_rev) * sc_unit_vec
                            '''

                            ###
                            # METHOD 2 for placing center points
                            ###

                            # Get the edge vertices, average them, then project onto the cable model
                            edge_vert_all_list = [item for sublist in bdry_e_list for item in sublist]
                            edge_vert_uniq_list = list(set(edge_vert_all_list))
                            edge_vert_co_list = [ob.data.vertices[item].co for item in edge_vert_uniq_list]
                            ave_co = Vector([0.0,0.0,0.0])
                            for item in edge_vert_co_list:
                                ave_co += item
                            ave_co /= len(edge_vert_co_list)
                            plane_ctr_pt = project_pt_line(this_sc.sc_pts[0],this_sc.sc_pts[1],ave_co)

                            # Make plane
                            pname = "sc_%02d_%02d_sg_%02d_B_sc_%02d_%02d_sg_%02d" % (pt_min,pt_max,i_sg+1,pt_min,pt_max,i_sg+2)
                            plane = MN_plane(name=pname, sides_ids=((pt_min,pt_max,i_sg+1),(pt_min,pt_max,i_sg+2)), edge_list=bdry_e_list, ctr_pt=plane_ctr_pt)

                            # Store the plane
                            this_sc.planes_sg_brdrs.append(plane)

                            # Store the number of segments
                            this_sc.n_seg = n_seg_rev

                else:

                    print("Number of segments: 1 (post revision)")

                    # Store the number of segments
                    this_sc.n_seg = 1

                    # Surface faces have already been stored for this section

            else:

                print("Number of segments: 1 (pre revision)")

                # Store the number of segments
                this_sc.n_seg = 1

                # Store the surface faces for this section
                surf_sg_f_list = [[]]
                for f in ob.data.polygons:
                    if f.select == True:
                        surf_sg_f_list[0].append(f.index)

                this_sc.surf_sg_f_list = surf_sg_f_list

            # Deselect all faces
            bpy.ops.object.mode_set(mode='EDIT') # Must be in edit mode
            bpy.ops.mesh.select_all(action='DESELECT')


    ###
    # Fix the names of the section planes
    ###

    # Make planes
    planes_done_list = []
    for sc_1 in mn_section_dict.values():
        for edge_list in sc_1.planes_sc_brdrs: 
            for plane in edge_list:
                if not plane.name in planes_done_list:

                    # Get the two sections on the plane's sides
                    ids_side_1 = sc_1.sc_id
                    side_ids = plane.sides_ids
                    if side_ids[0] == ids_side_1:
                        ids_side_2 = side_ids[1]
                    else:
                        ids_side_2 = side_ids[0]

                    # Get the other section
                    for sc in mn_section_dict.values():
                        if sc.sc_id == ids_side_2:
                            sc_2 = sc
                            break

                    # Which vertex are we at?
                    if ids_side_1[0] in ids_side_2:
                        vert_id = ids_side_1[0]
                    else:
                        vert_id = ids_side_1[1]

                    # For each section:
                    # if we are at the lower of the two vertices => segment 1
                    # if we are at the higher of the two => segment n_seg_#
                    if ids_side_1[0] == vert_id:
                        trip_1 = (ids_side_1[0],ids_side_1[1],1)
                    else:
                        trip_1 = (ids_side_1[0],ids_side_1[1],sc_1.n_seg)
                    if ids_side_2[0] == vert_id:
                        trip_2 = (ids_side_2[0],ids_side_2[1],1)
                    else:
                        trip_2 = (ids_side_2[0],ids_side_2[1],sc_2.n_seg)

                    # The new name of this plane
                    new_plane_name = "sc_%02d_%02d_sg_%02d" % min(trip_1,trip_2) + "_B_" + "sc_%02d_%02d_sg_%02d" % max(trip_1,trip_2)
                    print("Fixing name: " + str(plane.name) + " to: " + str(new_plane_name))
                    plane.name = new_plane_name

                    # Don't repeat this plane
                    planes_done_list.append(plane.name)


    ############################
    ############################
    # Make a surface plane object
    ############################
    ############################


    bpy.ops.object.mode_set(mode='OBJECT')
    # Select the original surface object
    # Deselect all
    for ob_tmp in bpy.data.objects:
        ob_tmp.select = False
    ob.select = True
    context.scene.objects.active = ob

    # Duplicate
    bpy.ops.object.duplicate()

    # Add the duplicate to MCell
    obj_surf = bpy.data.objects[ob.name + ".001"]
    obj_surf.name = obj_surf.name[:-4] + "_Surface"

    # Select the new object
    # Deselect all
    for ob_tmp in bpy.data.objects:
        ob_tmp.select = False
    obj_surf.select = True
    context.scene.objects.active = obj_surf

    # Add to MCell
    preserve_selection_use_operator(bpy.ops.mcell.model_objects_add, obj_surf)

    # Remove all existing MCell regions
    obj_surf.mcell.regions.remove_all_regions(context)

    # Update the object
    obj_surf.data.update()

    # Add each of the surfaces as a region
    planes_done_list = []
    for sc_id, sc in mn_section_dict.items():
        for i_sg in range(0,sc.n_seg):
            # Region name
            reg_name = "sc_%02d_%02d_sg_%02d_B_surf" % (sc_id[0],sc_id[1],i_sg+1)
    
            # Make region
            obj_surf.mcell.regions.add_region_by_name(context,reg_name)
            # Get region
            new_reg = obj_surf.mcell.regions.region_list[reg_name]

            # Assign faces
            new_reg.set_region_faces(obj_surf.data, sc.surf_sg_f_list[i_sg])

            # Prevent double counting
            planes_done_list.append(plane.name)

    # Update (but dont validate because who knows)
    obj_surf.data.update()


    ############################
    ############################
    # Make a segment plane object
    ############################
    ############################


    # Vertex, edge, face lists
    v_new_list = []
    v_id_dict = {}
    e_new_list = []
    f_new_list = []
    plane_f_dict = {}

    # Original vertex list
    ob_vert_list = [tuple(item.co) for item in ob.data.vertices]

    # Go through all the planes
    planes_done_list = []
    for sc_id, sc in mn_section_dict.items():
        # List of all planes
        all_plane_list = sc.planes_sg_brdrs
        for end_list in sc.planes_sc_brdrs:
            all_plane_list += end_list
        # Go through all planes
        for plane in all_plane_list:
            if not plane.name in planes_done_list:
                
                # Add this to the lists
                
                # List of unique vertex ids
                v_all_list = [item for sublist in plane.edge_list for item in sublist]
                v_unique_list = list(set(v_all_list))

                # Add to the vertex list
                for v in v_unique_list:
                    if not v in v_id_dict:
                        v_new_list.append(tuple(ob_vert_list[v]))
                        v_id_dict[v] = len(v_new_list)-1

                v_new_list.append(tuple(plane.ctr_pt))
                ctr_id = len(v_new_list) - 1

                # Add to the face list
                plane_f_dict[plane.name] = [] # Also, store what the face ids are in the new obj
                for e in plane.edge_list:
                    f_new_list.append([v_id_dict[e[0]],v_id_dict[e[1]],ctr_id])
                    plane_f_dict[plane.name].append(len(f_new_list)-1)

                # Add to the edge list
                e_new_list = []
                for e in plane.edge_list:
                    e_new_list.append([v_id_dict[e[0]],v_id_dict[e[1]]])
                for v in v_unique_list:
                    e_new_list.append([ctr_id,v_id_dict[v]])

                # Store as done
                planes_done_list.append(plane.name)

    # Make the object
    obj_name = ob.name + "_Segment"
    mesh_new = bpy.data.meshes.new(obj_name + "_mesh")
    mesh_new.from_pydata(v_new_list,e_new_list,f_new_list)
    # Validate and update
    mesh_new.validate(verbose=False) # Important! and i dont know why
    mesh_new.update()
    # Overwrite existing object
    obj_old = bpy.data.objects.get(obj_name)
    if obj_old:
        context.scene.objects.unlink(obj_old)
        bpy.data.objects.remove(obj_old)
    # New one
    obj_seg = bpy.data.objects.new(obj_name,mesh_new)
    # Link
    context.scene.objects.link(obj_seg)

    # Select the new object
    # Deselect all
    for ob_tmp in bpy.data.objects:
        ob_tmp.select = False
    obj_seg.select = True
    context.scene.objects.active = obj_seg

    # Add to MCell
    preserve_selection_use_operator(bpy.ops.mcell.model_objects_add, obj_seg)

    # Update the object
    obj_seg.data.update()

    # Add each of the surfaces as a region
    planes_done_list = []
    for sc_id, sc in mn_section_dict.items():
        # List of all planes
        all_plane_list = sc.planes_sg_brdrs
        for end_list in sc.planes_sc_brdrs:
            all_plane_list += end_list
        # Go through all planes
        for plane in all_plane_list:
            if not plane.name in planes_done_list:

                # Region name
                reg_name = plane.name
        
                # Make region
                obj_seg.mcell.regions.add_region_by_name(context,reg_name)
                # Get region
                new_reg = obj_seg.mcell.regions.region_list[reg_name]

                # Assign faces
                new_reg.set_region_faces(obj_seg.data, plane_f_dict[plane.name])

                # Prevent double counting
                planes_done_list.append(plane.name)

    # Update (but dont validate because who knows)
    obj_seg.data.update()

    print("> Finished: f_compartmentize_cyl")
    print("> Time: " + str(t_st[-1] - t_st[0]))

    return



