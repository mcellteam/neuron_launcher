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

# Function to write out to a file the number of segments in each section
def write_nseg_sec(fname):
    global mn_section_dict

    f = open(fname,'w')
    
    # Go through all sections
    ids = list(mn_section_dict.keys())
    n = len(ids)
    for i,sec_id in enumerate(ids):
        
        # Write the number of segments
        f.write(str(sec_id[0]) + " " + str(sec_id[1]) + " " + str(mn_section_dict[sec_id].n_seg))
        if i != n-1:
            f.write("\n")
    
    f.close()

# Function to write region names (surfaces and segments)
def write_nm_region_names_areas(fname):
    global nm_segment_dict

    f = open(fname,'w')
    
    # Write surface planes
    f.write('NM_Surface\n')
    n = len(list(nm_segment_dict.keys()))
    names_done = []
    for i_seg,seg_id in enumerate(list(nm_segment_dict.keys())):
        seg = nm_segment_dict[seg_id]
        if seg.plane_surf != -1: # segment has no plane facing the surface
            if not seg.plane_surf.name in names_done:
                f.write(seg.plane_surf.name + " " + str(seg.plane_surf.surf_area()) + "\n")
                names_done.append(seg.plane_surf.name)
        
    # Write segment planes
    f.write('NM_Segment\n')
    n = len(list(nm_segment_dict.keys()))
    names_done = []
    for i_seg,seg_id in enumerate(list(nm_segment_dict.keys())):
        seg = nm_segment_dict[seg_id]

        m = len(seg.planes_sg)
        for i_plane,plane_sg in enumerate(seg.planes_sg):
            if not plane_sg.name in names_done:
                
                f.write(plane_sg.name + " " + str(plane_sg.surf_area()))
            
                if not (i_seg == n-1 and i_plane == m-1):
                    f.write("\n")

                names_done.append(plane_sg.name)

    f.close()

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

        # Dictionary of segment index to face indexes of the surface
        self.surf_plane_f_list = {}

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

# Replace a mesh of data
def replace_mesh(obj_to_clear, rep_mesh):
    
    # Get the Bmesh
    if obj_to_clear.data.is_editmode:
        bm = bmesh.from_edit_mesh(obj_to_clear.data)
    else:
        bm = bmesh.new()
        bm.from_mesh(obj_to_clear.data)
    
    # Clear everything!
    bm.clear()

    # Replace
    bm.from_mesh(rep_mesh)

    # Update object
    if bm.is_wrapped:
        bmesh.update_edit_mesh(obj_to_clear.data)
    else:
        bm.to_mesh(obj_to_clear.data)
        obj_to_clear.data.update()

# Clear a mesh of all data - HOW THE FUCK DO YOU DO THIS WITHOUT BMESH????
def clear_mesh(obj_to_clear):

    # Get the Bmesh
    if obj_to_clear.data.is_editmode:
        bm = bmesh.from_edit_mesh(obj_to_clear.data)
    else:
        bm = bmesh.new()
        bm.from_mesh(obj_to_clear.data)

    # Clear everything!
    bm.clear()

    # Update object
    if bm.is_wrapped:
        bmesh.update_edit_mesh(obj_to_clear.data)
    else:
        bm.to_mesh(obj_to_clear.data)
        obj_to_clear.data.update()

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

if __name__ == "__main__":

    print("> Running: f_compartmentize_fast")

    # Get the working dir from the filepath
    global DIR
    DIR = "/Users/oernst/Research/cnl/neuron_mcell/grant_proposal/figures_blender/"
    #swc_filename = os.path.basename(swc_filepath)
    #DIR = swc_filepath[:-len(swc_filename)]
    #print("Working directory: " + str(DIR))

    # Dictionary of sections
    global mn_section_dict
    mn_section_dict = {}

    # Time
    t_st = []
    t_st.append(time.time())

    # Get data from the SWC file
    get_connections(DIR + "cable_model_cut.swc")

    # Get the active object
    ob_list = bpy.context.selected_objects

    if len(ob_list) != 1:
        sys.exit()
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
            sec.select_region_faces(bpy.context)

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
    bpy.context.scene.objects.active = ob

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

    # Number of segments per length
    n_seg_plen = 1.5

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

            if n_seg != 1:
                    
                print("Making list of dividing points")

                # Make a list of all the dividing points on the cable section
                d_sg_length = sc_length / n_seg
                sc_dividing_pt_list = []
                for i in range(1,n_seg):
                    sc_dividing_pt_list.append(this_sc.sc_pts[0] + i * d_sg_length * sc_unit_vec)

                # List of dividing lengths from sc_pts[0]
                sc_dividing_len_list = [0.0]
                for pt in sc_dividing_pt_list:
                    sc_dividing_len_list.append((pt-this_sc.sc_pts[0]).length)
                sc_dividing_len_list.append(2.0*sc_length)

                print("Dividing lengths")
                print(sc_dividing_len_list)

                # List of face id's in each segment
                sg_f_list = []
                for i in range(0,n_seg):
                    sg_f_list.append([])

                # Go through all faces in this section; project them to the cable section
                sec.select_region_faces(bpy.context)
                # bpy.ops.object.mode_set(mode='OBJECT')
                for f in ob.data.polygons:
                    if f.select == True:
                        proj_pt = project_pt_line(this_sc.sc_pts[0],this_sc.sc_pts[1],Vector(f.center))
                        proj_dist = (proj_pt - this_sc.sc_pts[0]).length
                        i_sg = bisect.bisect(sc_dividing_len_list, proj_dist) - 1

                        # Store face
                        sg_f_list[i_sg].append(f.index)

                # Go through each segment, get the boundary edge
                sg_bdry_e_list = []
                for i in range(0,n_seg):
                    sg_bdry_e_list.append([])
                for i_sg, f_list in enumerate(sg_f_list):
                    # Deselect all
                    # bpy.ops.object.mode_set(mode='EDIT')
                    bpy.ops.mesh.select_all(action='DESELECT')

                    # Select all faces
                    for f in f_list:
                        ob.data.polygons[f].select = True

                    # Get the border loop
                    bpy.ops.mesh.region_to_loop()
                    for e in ob.data.edges:
                        if e.select == True:
                            sg_bdry_e_list[i_sg].append(e.index)

                # Determine boundaries
                for i_sg in range(0,n_seg-1):
                    bdry_list_1 = sg_bdry_e_list[i_sg]
                    bdry_list_2 = sg_bdry_e_list[i_sg+1]

                    # Check for overlap
                    bdry_overlap = list(set(bdry_list_1).intersection(set(bdry_list_2)))

                    if len(bdry_overlap) > 0:

                        # Convert to edge pair list
                        bdry_e_list_tmp = [ob.data.edges[item].vertices for item in bdry_overlap]
                        bdry_e_list = [[int(ids[0]),int(ids[1])] for ids in bdry_e_list_tmp]

                        # Make plane
                        pname = "sc_%02d_%02d_sg_%02d_B_sc_%02d_%02d_sg_%02d" % (pt_min,pt_max,i_sg+1,pt_min,pt_max,i_sg+2)
                        plane = MN_plane(name=pname, sides_ids=((pt_min,pt_max,i_sg+1),(pt_min,pt_max,i_sg+2)), edge_list=bdry_e_list, ctr_pt=sc_dividing_pt_list[i_sg])

                        # Store the plane
                        this_sc.planes_sg_brdrs.append(plane)

            else:

                # Only one segment for the whole section
                sec.select_region_faces(bpy.context)

                # Assign surface faces
                sg_f_list = [[]]
                for f in ob.data.polygons:
                    if f.select == True:
                        sg_f_list[0].append(f.index)

            # bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='DESELECT')

    # Make planes
    ob_vert_list = [tuple(item.co) for item in ob.data.vertices]
    planes_done_list = []
    for sc in mn_section_dict.values():
        for end_list in sc.planes_sc_brdrs: 
            for plane in end_list:
                if not plane.name in planes_done_list:
                    plane.make_separate_object(ob_vert_list)
                    planes_done_list.append(plane.name)

    '''
    ###
    # Mesh errr thang
    ###

    # Global list of tet face marker to plane it belongs to
    global face_marker_plane_dict
    face_marker_plane_dict = {}

    mesh_built = make_tetgen_mesh()

    # Make it an object
    #print("> Creating a Blender object from a TetGen mesh....")
    #new_obj_meshpy("tetgen", mesh_built);
    # mesh_built.save_elements(DIR+"tetgen_elements")

    # Get the coordinates of the tet points
    mesh_vert_co_list = [tuple(item) for item in list(mesh_built.points)]

    # List of tets to vertex indices
    mesh_tet_vert_list = list(mesh_built.elements)

    # List of neighbours of each tet
    mesh_tet_nghbr_list = list(mesh_built.neighbors)

    # List of face vertex ids to face markers
    mesh_face_vert_marker_dict = dict(mesh_built.face_vertex_indices_to_face_marker)

    # Tet attributes - every index is a unique (but unknown :'( ) region
    mesh_tet_att = list(mesh_built.element_attributes)

    # # # Reclaim memory
    del mesh_built

    ## DEBUG = UNIMPORTANT TO FUNCTION
    # All the unique attributes present
    unique_tet_att = list(set(mesh_tet_att))
    print("> > Unique tet attributes present: ")
    print(unique_tet_att)
    print(len(unique_tet_att))
    ##

    # Dictionary of section attribute to tet id
    sc_att_tet_dict = {}
    for i_tet, sc_att in enumerate(mesh_tet_att):
        if int(sc_att) in sc_att_tet_dict:
            sc_att_tet_dict[int(sc_att)].append(i_tet)
        else:
            sc_att_tet_dict[int(sc_att)] = [i_tet]

    # Figure out what attribute indices correspond to which section
    sc_tet_dict = {}
    for sc_att in sc_att_tet_dict.keys():
        sc_tets = sc_att_tet_dict[sc_att]
        
        # Find surface tets
        DONE = False
        for tet in sc_tets:
            if -1 in mesh_tet_nghbr_list[tet]:
                # Get its verts
                tet_verts = mesh_tet_vert_list[tet]
                # Get the triplets of face vertex indices
                for triplet in list(itertools.combinations(tet_verts,3)):
                    # Get the face marker
                    face_marker = mesh_face_vert_marker_dict[frozenset(triplet)]
                    # Is it a surface?
                    if face_marker < 0:
                        # What plane corresponds to this marker
                        plane = face_marker_plane_dict[face_marker]
                        # What is the section key for this plane
                        sc_id = plane.sides_ids[0]
                        # Store
                        sc_tet_dict[sc_id] = sc_tets
                        # Stop for this section
                        DONE = True
                        break
                            
                if DONE == True:
                    break

    ## DEBUG = UNIMPORTANT TO FUNCTION
    # All attributes that were assigned
    print("> > All att that were assigned were assigned to these sections:")
    print(sc_tet_dict.keys())
    print(len(sc_tet_dict.keys()))
    ##

    # Time
    t_st.append(time.time())
    print("> Created Tetgen Mesh: time: " + str(t_st[-1]-t_st[-2]))

    ###
    # Segment each of the sections
    ###

    # Number of segments per length
    n_seg_plen = 1.5

    for sec_id in list(mn_section_dict.keys()):
        
        # if sec_id in sc_tet_dict: # Why should this not be the case? This occurs for some reason....
        
        sec = mn_section_dict[sec_id]

        # Make the dividing planes
        segment_meshpy(mesh_vert_co_list, mesh_tet_vert_list, sc_tet_dict[sec_id], mesh_tet_nghbr_list, mesh_face_vert_marker_dict, sec, n_seg_plen)

    # Time
    t_st.append(time.time())
    print("> Created segment borders for each section: time: " + str(t_st[-1]-t_st[-2]))

    # # # Free memory
    del mesh_vert_co_list
    del mesh_tet_vert_list
    del mesh_tet_nghbr_list
    del mesh_face_vert_marker_dict
    del mesh_tet_att

    ###
    # After making all the segments: Convert segment->section boundaries into segment->segment boundaries!
    ###

    # Time
    t_st.append(time.time())

    # Go through all of the segments
    for nm_seg_id in nm_segment_dict.keys():
        nm_seg = nm_segment_dict[nm_seg_id]

        # Go through all of the section bounding planes
        for plane_sc in nm_seg.planes_sc:

            # What's on the other side of this plane?
            sides_ids = plane_sc.sides_ids
            if len(sides_ids[0]) == 2 and len(sides_ids[1]) == 3:
                other_side_sc = sides_ids[0]
            elif len(sides_ids[1]) == 2 and len(sides_ids[0]) == 3:
                other_side_sc = sides_ids[1]
            else:
                # Not a section -> segment plane - somehow?
                print("Warning! This isn't a segment->section plane! Should this be possible?")
            
            # Get all the segments that are in the other side's section
            for other_nm_seg_id in nm_segment_dict.keys():

                # Obviously, check yourself before you wreck yourself
                if nm_seg_id != other_nm_seg_id:
                    # Is this segment in the section we care about
                    if other_nm_seg_id[0:2] == other_side_sc:

                        # Get the segment
                        other_nm_seg = nm_segment_dict[other_nm_seg_id]

                        # Check all the planes in this segment for overlap
                        for other_plane_sc in other_nm_seg.planes_sc:
                            
                            # What's on the other side of this plane?
                            other_sides_ids = other_plane_sc.sides_ids
                            if len(other_sides_ids[0]) == 2 and len(other_sides_ids[1]) == 3:
                                other_other_side_sc = other_sides_ids[0]
                            elif len(other_sides_ids[1]) == 2 and len(other_sides_ids[0]) == 3:
                                other_other_side_sc = other_sides_ids[1]
                            else:
                                # Not a section -> segment plane - somehow?
                                print("Warning! This isn't a segment->section plane! Should this be possible?")

                            # Is the section on the other side of this plane to original one?
                            if nm_seg_id[0:2] == other_other_side_sc:
                                
                                # The would be plane id
                                min_id = min(nm_seg_id,other_nm_seg_id)
                                max_id = max(nm_seg_id,other_nm_seg_id)
                                p_sides_ids = (min_id,max_id)

                                # Make sure we didn't check this already
                                if not p_sides_ids in [plane0.sides_ids for seg0 in nm_segment_dict.values() for plane0 in seg0.planes_sg]:
                                
                                    # Check for overlap
                                    vert_list, face_list = plane_sc.overlap(other_plane_sc)
                                    if vert_list != None:
                                        
                                        # Yes! there's overlap
                                        # Make yet another plane
                                        p_name = 'sc_%02d_%02d_sg_%02d'%min_id + '_B_' + 'sc_%02d_%02d_sg_%02d'%max_id
                                        p_sides_names = ('sc_%02d_%02d_sg_%02d'%min_id, 'sc_%02d_%02d_sg_%02d'%max_id)
                                        plane = MN_plane(name=p_name, sides_names=p_sides_names, sides_ids=p_sides_ids, vert_list=vert_list, face_list=face_list, CONV=True)
                                        
                                        # Add the plane to both segments
                                        nm_seg.planes_sg.append(plane)
                                        other_nm_seg.planes_sg.append(plane)


    # Time
    t_st.append(time.time())
    print("> Fixed correspondence of section-section planes: time: " + str(t_st[-1]-t_st[-2]))

    ###
    # Create a single object out of all the plane surfaces with MCell regions
    ###
    
    plane_surf_vert_stack = []
    plane_surf_face_stack = []
    
    # Go through all segments
    plane_sides_ids_done = [] # Make sure not to double count any planes
    for nm_seg_id in nm_segment_dict.keys():
        nm_seg = nm_segment_dict[nm_seg_id]
        
        if nm_seg.plane_surf != -1: # No surface for this segment (definitely possible!)
            
            if not nm_seg.plane_surf.sides_ids in plane_sides_ids_done:
                
                plane_surf_vert_stack += [nm_seg.plane_surf.vert_list]
                plane_surf_face_stack += [nm_seg.plane_surf.face_list]

                # Prevent double counting
                plane_sides_ids_done.append(nm_seg.plane_surf.sides_ids)

    # Correct the lists
    plane_surf_vert_list, plane_surf_face_list = stack_lists(plane_surf_vert_stack, plane_surf_face_stack)

    # Make the plane
    nm_surface_plane = MN_plane(name="NM_Surface",vert_list = plane_surf_vert_list, face_list = plane_surf_face_list, CONV=True)

    # Make make the plane
    obj_nm_surface_plane = nm_surface_plane.make_plane()

    # Set active
    context.scene.objects.active = obj_nm_surface_plane
    obj_nm_surface_plane.select = True

    # Add to MCell
    preserve_selection_use_operator(bpy.ops.mcell.model_objects_add, obj_nm_surface_plane)

    # Add each of the segments as a region
    face_idx_ctr = 0
    plane_sides_ids_done = [] # Make sure not to double count any planes
    for nm_seg_id in nm_segment_dict.keys():
        nm_seg = nm_segment_dict[nm_seg_id]
        
        if nm_seg.plane_surf != -1: # No surface for this segment (definitely possible!)
            
            if not nm_seg.plane_surf.sides_ids in plane_sides_ids_done:
                
                # Region name
                reg_name = nm_seg.plane_surf.name

                # Make region
                obj_nm_surface_plane.mcell.regions.add_region_by_name(context,reg_name)
                # Get region
                new_reg = obj_nm_surface_plane.mcell.regions.region_list[reg_name]
                
                # Assign faces
                face_ids = range(face_idx_ctr,face_idx_ctr+len(nm_seg.plane_surf.face_list))
                face_idx_ctr += len(nm_seg.plane_surf.face_list)
                new_reg.set_region_faces(obj_nm_surface_plane.data, face_ids)

                # Prevent double counting
                plane_sides_ids_done.append(nm_seg.plane_surf.sides_ids)

    # Update (but dont validate because who knows)
    obj_nm_surface_plane.data.update()

    # Deselect
    obj_nm_surface_plane.select = False

    # Time
    t_st.append(time.time())
    print("> Created surface plane: time: " + str(t_st[-1]-t_st[-2]))

    ###
    # Create a single object out of all the segments surfaces with MCell regions
    ###

    seg_surf_vert_stack = []
    seg_surf_face_stack = []

    # Go through all segments
    plane_sides_ids_done = [] # Make sure not to double count any planes
    for nm_seg_id in nm_segment_dict.keys():
        nm_seg = nm_segment_dict[nm_seg_id]
        
        # Go through all it's segment bounding planes
        for plane_sg in nm_seg.planes_sg:
            
            if not plane_sg.sides_ids in plane_sides_ids_done:
                seg_surf_vert_stack += [plane_sg.vert_list]
                seg_surf_face_stack += [plane_sg.face_list]
                
                # Prevent double counting
                plane_sides_ids_done.append(plane_sg.sides_ids)

    # Correct the lists
    seg_surf_vert_list, seg_surf_face_list = stack_lists(seg_surf_vert_stack, seg_surf_face_stack)

    '''

    '''
    # TEMP: throw out pieces outside the desired range
    tmp_vert_list = []
    tmp_dict = {}
    tmp_face_list = []
    # What verts are allowed
    for i_v,v in enumerate(seg_surf_vert_list):
        if v[1] > -9421.5459-1.0 and v[1] < -9421.5459+1.0 and v[2] > 315.18835-1.0 and v[2] < 315.18835+1.0:
            tmp_vert_list.append(v)
            tmp_dict[i_v] = len(tmp_vert_list) - 1
    # What faces are allowed
    tmp_keys = list(tmp_dict.keys())
    for f in seg_surf_face_list:
        if f[0] in tmp_keys and f[1] in tmp_keys and f[2] in tmp_keys:
            tmp_face_list.append([tmp_dict[f[0]],tmp_dict[f[1]],tmp_dict[f[2]]])

    print("> Constructed temp vert and face list:")
    print(len(tmp_vert_list))
    print(len(tmp_face_list))

    nm_seg_plane = MN_plane(name="NM_Segment",vert_list = tmp_vert_list, face_list = tmp_face_list, CONV=True)
    '''

    '''

    # Make the plane
    nm_seg_plane = MN_plane(name="NM_Segment",vert_list = seg_surf_vert_list, face_list = seg_surf_face_list, CONV=True)

    # Make make the plane
    obj_nm_seg_plane = nm_seg_plane.make_plane()

    # Set active
    context.scene.objects.active = obj_nm_seg_plane
    obj_nm_seg_plane.select = True

    # Add to MCell
    preserve_selection_use_operator(bpy.ops.mcell.model_objects_add, obj_nm_seg_plane)

    # Add each of the segments as a region
    face_idx_ctr = 0
    plane_sides_ids_done = [] # Make sure not to double count any planes
    for nm_seg_id in nm_segment_dict.keys():
        nm_seg = nm_segment_dict[nm_seg_id]
        
        # Go through all it's segment bounding planes
        for plane_sg in nm_seg.planes_sg:

            if not plane_sg.sides_ids in plane_sides_ids_done:
                
                # Region name
                reg_name = plane_sg.name
                
                # Make region
                obj_nm_seg_plane.mcell.regions.add_region_by_name(context,reg_name)
                # Get region
                new_reg = obj_nm_seg_plane.mcell.regions.region_list[reg_name]

                # Assign faces
                face_ids = range(face_idx_ctr,face_idx_ctr+len(plane_sg.face_list))
                face_idx_ctr += len(plane_sg.face_list)
                new_reg.set_region_faces(obj_nm_seg_plane.data, face_ids)

                # Prevent double counting
                plane_sides_ids_done.append(plane_sg.sides_ids)

    # Update (but dont validate because who knows)
    obj_nm_seg_plane.data.update()

    # Deselect
    obj_nm_seg_plane.select = False

    # Time
    t_st.append(time.time())
    print("> Created segment plane: time: " + str(t_st[-1]-t_st[-2]))


    ###
    # Fix normals
    ###

    ###
    # Surface
    ###

    # Set active
    context.scene.objects.active = obj_nm_surface_plane
    obj_nm_surface_plane.select = True

    # Edit mode
    bpy.ops.object.mode_set(mode='EDIT')

    # Go through all the regions in the object
    for reg in obj_nm_surface_plane.mcell.regions.region_list:
        # Select region
        reg.select_region_faces(context)

        # Fix normals
        bpy.ops.mesh.normals_make_consistent(inside=False)

        # Deselect all
        bpy.ops.mesh.select_all(action='DESELECT')

    # Object mode
    bpy.ops.object.mode_set(mode='OBJECT')

    # Deselect
    obj_nm_surface_plane.select = False

    ###
    # Segments
    ###

    # Set active
    context.scene.objects.active = obj_nm_seg_plane
    obj_nm_seg_plane.select = True
    
    # Edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    
    # Go through all the regions in the object
    for reg in obj_nm_seg_plane.mcell.regions.region_list:
        # Select region
        reg.select_region_faces(context)
        
        # Fix normals
        bpy.ops.mesh.normals_make_consistent(inside=False)
        
        # Check that we are pointing toward the bigger index
        reg_name = reg.name
        id_min_min = int(reg_name[3:5])
        id_min_max = int(reg_name[6:8])
        id_max_min = int(reg_name[20:22])
        id_max_max = int(reg_name[23:25])
        if id_min_min == id_max_min and id_min_max == id_max_max:
            # Segment boundary
            sc_pts = mn_section_dict[(id_min_min,id_min_max)].sc_pts
            vec_check = sc_pts[1]-sc_pts[0]
        else:
            # Section boundary
            all_ids = [id_min_min, id_min_max, id_max_min, id_max_max]
            shared_id = list(set([item for item in all_ids if all_ids.count(item) == 2]))[0]
            uniq_ids = [item for item in all_ids if all_ids.count(item) == 1]
            max_uniq_id = max(uniq_ids)
            min_id = min(shared_id,max_uniq_id)
            max_id = max(shared_id,max_uniq_id)
            sc_pts = mn_section_dict[(min_id,max_id)].sc_pts
            vec_check = sc_pts[1]-sc_pts[0]

        # Finally, check any face - dot product should be positive!
        for face in obj_nm_seg_plane.data.polygons:
            if face.select == True:
                fSel = face.vertices
                break

        v_fSel = [obj_nm_seg_plane.data.vertices[v].co for v in fSel]
        vec_fSel = (v_fSel[1]-v_fSel[0]).cross(v_fSel[2]-v_fSel[1])
        # Dot product
        dp = vec_fSel.dot(vec_check)
        if dp < 0:
            # Flip?!
            bpy.ops.mesh.flip_normals()

        # Deselect all
        bpy.ops.mesh.select_all(action='DESELECT')

    # Object mode
    bpy.ops.object.mode_set(mode='OBJECT')

    # Deselect
    obj_nm_surface_plane.select = False

    ###
    # Write out some relevant information
    ###

    write_nseg_sec(DIR+"nseg_sec.txt")

    write_nm_region_names_areas(DIR+"region_names_areas.txt")

    ###
    # Make individual objects out of all of the planes
    ###
    '''

    '''
        
    # Go through all segments
    plane_sides_ids_done = [] # Make sure not to double count any planes
    for nm_seg_id in nm_segment_dict.keys():
        nm_seg = nm_segment_dict[nm_seg_id]

        # Go through all it's segment bounding planes
        for plane_sg in (nm_seg.planes_sg + [nm_seg.plane_surf]):
            
            if plane_sg != -1: # No surface for this compartment
                
                if not plane_sg.sides_ids in plane_sides_ids_done:

                    plane_sg.make_plane()

                    # Counter
                    plane_sides_ids_done.append(plane_sg.sides_ids)

    # Time
    t_st.append(time.time())
    print("> Created objects: time: " + str(t_st[-1]-t_st[-2]))
    
    '''



    print("> Finished: f_compartmentize")
    print("> Time: " + str(t_st[-1] - t_st[0]))








