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

# MeshPy
from meshpy.tet import MeshInfo, build, Options

# Time
import time

# Function to write out to a file the number of segments in each section
def write_nseg_sec(fname):
    global nm_section_dict

    f = open(fname,'w')
    
    # Go through all sections
    ids = list(nm_section_dict.keys())
    n = len(ids)
    for i,sec_id in enumerate(ids):
        
        # Write the number of segments
        f.write(str(sec_id[0]) + " " + str(sec_id[1]) + " " + str(nm_section_dict[sec_id].n_seg))
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

# Function to make a tetgen mesh from the global nm_section list
def make_tetgen_mesh():
    global nm_section_dict
    global face_marker_plane_dict

    # Make single volume
    vert_stack = []
    face_stack = []
    face_marker_stack = []
    face_marker_brdr = 1 # are positive
    face_marker_surf = -1 # are negative

    # Make sure we don't double count any planes
    # Sides ids to face marker
    plane_sides_ids_done = {}

    # Go through all the sections
    for sec_id in nm_section_dict.keys():
        sec = nm_section_dict[sec_id]

        # Check that there are any borders
        if len(sec.planes_sc_brdrs[0]) == 0 and len(sec.planes_sc_brdrs[1]) == 0 and sec.plane_surf == -1:
            print("Error! No bounding planes for this section!")

        # Add the surface plane...
        if not sec.plane_surf.sides_ids in plane_sides_ids_done:
            vert_stack += [sec.plane_surf.vert_list]
            face_stack += [sec.plane_surf.face_list]
            face_marker_stack += [len(sec.plane_surf.face_list)*[face_marker_surf]]
            face_marker_plane_dict[face_marker_surf] = sec.plane_surf
            sec.face_marker_dict[face_marker_surf] = -1

            plane_sides_ids_done[sec.plane_surf.sides_ids] = face_marker_surf
        
            face_marker_surf -= 1
        else:
            # Do this for EVERY section
            sec.face_marker_dict[plane_sides_ids_done[sec.plane_surf.sides_ids]] = -1

        # And the borders...
        for i,brdr in enumerate(sec.planes_sc_brdrs):
            for j,plane in enumerate(brdr):
                if not plane.sides_ids in plane_sides_ids_done:
                    vert_stack += [plane.vert_list]
                    face_stack += [plane.face_list]
                    # Face marker
                    face_marker_stack += [len(plane.face_list)*[face_marker_brdr]]
                    face_marker_plane_dict[face_marker_brdr] = plane
                    sec.face_marker_dict[face_marker_brdr] = (i,j)
                    
                    '''
                    if nm_sec.sc_id == (1,3):
                        global shared_3_27
                        PLANE_SHARED = True
                        for v in shared_3_27:
                            if not v in plane.vert_list:
                                PLANE_SHARED = False
                                break
                    
                        if PLANE_SHARED:
                            print("THIS IS THE SHARED PLANE: marker: " + str(face_marker_brdr))
                        
                        if triplet[0] in shared_3_27 and triplet[1] in shared_3_27 and triplet[2] in shared_3_27:
                            print("Trip: " + str(triplet) + " is shared; has face marker: " + str(face_marker))
                    '''
                    
                    plane_sides_ids_done[plane.sides_ids] = face_marker_brdr
                
                    face_marker_brdr += 1
                else:
                    # Do this for EVERY section
                    sec.face_marker_dict[plane_sides_ids_done[plane.sides_ids]] = (i,j)

    # Create a correctly indexed closed volume
    vert_list, face_list, face_marker_list = stack_lists(vert_stack,face_stack,face_marker_stack=face_marker_stack)
    
    # Make the tetgen mesh
    mesh_info = MeshInfo()
    
    # Points
    mesh_info.set_points(vert_list)
    
    # Faces
    mesh_info.set_facets(face_list, markers = face_marker_list)
    
    # --- TEMP ---
    '''
    # Make an object from the surface we are planning to tetrahedronalize
    mesh_new = bpy.data.meshes.new("pre_tet_mesh")
    mesh_new.from_pydata(vert_list,[],face_list)
    mesh_new.validate(verbose=False) # Important! and i dont know why
    mesh_new.update()
    obj_new = bpy.data.objects.new("pre_tet",mesh_new)
    context.scene.objects.link(obj_new)
    # return
    '''
    # --- FIN TEMP ---

    # Tetrahedralize
    # Options:
    # neighout = Write out neighbors
    # facesout = Write out faces
    # edgesout = Write out edges
    # regionattrib = Write out element_attributes = unique id for every tet in a distinct volume
    # nobisect = Dont alter surface
    print("> Starting TetGen")
    opts = Options(switches='pq', neighout = True, facesout = True, edgesout = True, regionattrib = True, verbose = True, docheck = True)
    mesh_built = build(mesh_info, options=opts)
    print("> Finished TetGen successfully")

    return mesh_built



# Function to take vertex list, face list and create a single vert/face/edge list
def stack_lists(vert_stack, face_stack, edge_stack=None, MAKE_EDGE_LIST=False, face_marker_stack=None, RETURN_DICT=False):
    vert_list = []
    face_list = []
    if MAKE_EDGE_LIST == True:
        edge_list = []
    
    if face_marker_stack != None:
        # Faces are not created/eliminated NOR re-indexed - so just flatten the list!
        face_marker_list = [item for sublist in face_marker_stack for item in sublist]

    if RETURN_DICT == True:
        # Return a dictionary of indices for each vertex list in the stack to the new vertex list
        idx_dict_ret = []
        for i in range(0,len(vert_stack)):
            idx_dict_ret.append({})
        
    # Go through all the vertex lists
    for i_stack,vert_list_in in enumerate(vert_stack):
        
        # Go through all verts
        idx_dict = {}
        for i_v,v in enumerate(vert_list_in):
            # Does vertex already exist? (Duplicate vertex)
            if v in vert_list:
                idx_dict[i_v] = vert_list.index(v)
            else:
                # Append this vertex
                idx_dict[i_v] = len(vert_list)
                if RETURN_DICT == True:
                    idx_dict_ret[i_stack][i_v] = len(vert_list)

                vert_list.append(v)

        # Face list / marker list
        face_list_in = face_stack[i_stack]

        # Re-index the faces
        face_list_in_ri = [[idx_dict[v] for v in f] for f in face_list_in]

        # Append the faces
        face_list += face_list_in_ri

        # Edge list
        if MAKE_EDGE_LIST == True:
            # Are edges provided?
            if edge_stack != None:
                # Get the edge list
                edge_list_in = edge_stack[i_stack]
                
                # Re-index the edges
                edge_list_in_ri = [[idx_dict[v] for v in e] for e in edge_list_in]

                # Append the edges
                edge_list += edge_list_in_ri
            # Create an edge list from the faces
            else:
                edge_list_make = []
                for f in face_list_in_ri:
                    edge_list_make.append([f[0],f[1]])
                    edge_list_make.append([f[0],f[2]])
                    edge_list_make.append([f[1],f[2]])

                # Append
                edge_list += edge_list_make

    if MAKE_EDGE_LIST == True:
        if face_marker_stack != None:
            if RETURN_DICT == True:
                return vert_list, face_list, face_marker_list, edge_list, idx_dict_ret
            else:
                return vert_list, face_list, face_marker_list, edge_list
        else:
            if RETURN_DICT == True:
                return vert_list, face_list, edge_list, idx_dict_ret
            else:
                return vert_list, face_list, edge_list
    else:
        if face_marker_stack != None:
            if RETURN_DICT == True:
                return vert_list, face_list, face_marker_list, idx_dict_ret
            else:
                return vert_list, face_list, face_marker_list
        else:
            if RETURN_DICT == True:
                return vert_list, face_list, idx_dict_ret
            else:
                return vert_list, face_list

# Class for a section
class NM_section:
    
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

        # Surface plane
        self.plane_surf = -1

        # Dictionary for tetgen face markers to indeces of planes
        self.face_marker_dict = {}

# Class for a segment
class NM_segment:

    # Init
    def __init__(self, *args, **kwargs):
        # Segment name
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = ""

        # Tuple of three (two sc + 1 sg) ids that define this segment
        if 'sg_id' in kwargs:
            self.sg_id = kwargs['sg_id']
        else:
            self.sg_id = (-1,-1,-1)

        # List of segment bounding planes
        self.planes_sg = []
        
        # List of section bounding planes
        self.planes_sc = []

        # Surface plane
        self.plane_surf = -1

# Class for a plane
class NM_plane:
    
    # Compute the surface area
    def surf_area(self):
        sa = 0.0
        for f in self.face_list:
            sa += 0.5*( (Vector(self.vert_list[f[1]]) - Vector(self.vert_list[f[0]])).cross(Vector(self.vert_list[f[2]]) - Vector(self.vert_list[f[0]])) ).length
        
        return sa
    
    # Make it into a blender object
    def make_plane(self):
        print("> > Making: " + str(self.name) + " into a Blender object.")
        print("> > Number of verts: " + str(len(self.vert_list)))
        print("> > Number of edges: " + str(len(self.edge_list)))
        print("> > Number of faces: " + str(len(self.face_list)))
        
        # Make an obj
        mesh_new = bpy.data.meshes.new(self.name + "_mesh")
        
        mesh_new.from_pydata(self.vert_list, self.edge_list, self.face_list)
        
        # Validate and update
        mesh_new.validate(verbose=False) # Important! and i dont know why
        mesh_new.update()
        
        # Delete an old object if it exists
        '''
        if bpy.data.objects.get(self.name) is not None:
            obj_del = bpy.data.objects.get(self.name)
            
            # Delete MCell obj if it exists
            mcell_obj_list = context.scene.mcell.model_objects.object_list
            mcell_obj_names = [item.name for item in mcell_obj_list]
            if self.name in mcell_obj_names:
            
                # Delete the object
                preserve_selection_use_operator(bpy.ops.mcell.model_objects_remove, obj_del)

            obj_del.select = True
            bpy.ops.object.delete()
        '''
        
        '''
        # Don't delete - instead, if the object exists, just clear the mesh
        if bpy.data.objects.get(self.name) is not None:
            obj_new = bpy.data.objects.get(self.name)
            replace_mesh(obj_new, mesh_new)
        
        else:
        '''
        
        # Overwrite existing object
        obj_old = bpy.data.objects.get(self.name)
        if obj_old:
            bpy.context.scene.objects.unlink(obj_old)
            bpy.data.objects.remove(obj_old)
            
        # New one
        obj_new = bpy.data.objects.new(self.name,mesh_new)
        
        print("> > Number of verts in new obj: " + str(len(obj_new.data.vertices)))
        print("> > Number of edges in new obj: " + str(len(obj_new.data.edges)))
        print("> > Number of faces in new obj: " + str(len(obj_new.data.polygons)))

        # Link
        bpy.context.scene.objects.link(obj_new)

        return obj_new

    # Check for overlap with another plane
    def overlap(self, op):
        # Convert all of my faces into sorted point triplets
        self_face_trip = []
        for f in self.face_list:
            # Figure out the triplet of points
            pt_trip = [self.vert_list[f[0]],self.vert_list[f[1]],self.vert_list[f[2]]]
            pt_trip.sort()
            self_face_trip.append(pt_trip)
        
        # And the other guys
        op_face_trip = []
        for f in op.face_list:
            # Figure out the triplet of points
            pt_trip = [op.vert_list[f[0]],op.vert_list[f[1]],op.vert_list[f[2]]]
            pt_trip.sort()
            op_face_trip.append(pt_trip)
        
        # Go through all faces in one, search
        shared_face_trip = []
        for trip in self_face_trip:
            # Easy check
            if trip in op_face_trip:
                shared_face_trip.append(trip)
        # Hard check
        '''
        for op_trip in op_face_trip:
            d0 = (Vector(trip[0])-Vector(op_trip[0])).length
            d1 = (Vector(trip[1])-Vector(op_trip[1])).length
            d2 = (Vector(trip[2])-Vector(op_trip[2])).length
            if max(d0,d1,d2) < 1e-5:
                shared_face_trip.append(trip)
        '''
    
        if len(shared_face_trip) == 0:
            return None, None
        
        # Convert to vert list / face list
        v_list = list(set([item for sublist in shared_face_trip for item in sublist]))
        f_list = []
        for trip in shared_face_trip:
            f_list.append([v_list.index(trip[0]),v_list.index(trip[1]),v_list.index(trip[2])])
                
        return v_list, f_list
    
    # Convert full list of vert, faces to this one's verts, edges, faces
    def conv_lists(self, verts0, faces0):
        # Clear
        self.vert_list = []
        self.edge_list = []
        self.face_list = []
        # Dictionary of old index to new
        dict_o_i = {}
        for f in faces0:
            # Vert + dict
            for v in f:
                if not v in dict_o_i:
                    dict_o_i[v] = len(self.vert_list)
                    self.vert_list.append(tuple(verts0[v]))
            # Edges
            self.edge_list.append([dict_o_i[f[0]],dict_o_i[f[1]]])
            self.edge_list.append([dict_o_i[f[0]],dict_o_i[f[2]]])
            self.edge_list.append([dict_o_i[f[1]],dict_o_i[f[2]]])
            
            # Face
            self.face_list.append([dict_o_i[f[0]],dict_o_i[f[1]],dict_o_i[f[2]]])

    # Init
    def __init__(self, *args, **kwargs):
        # Name
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = ""
        
        # Names for what's on either side
        if 'sides_names' in kwargs:
            self.sides_names = kwargs['sides_names']
        else:
            self.sides_names = ("","")
        
        # Triplets (two sc + 1 sg) ids for what's on either side, lower first
        if 'sides_ids' in kwargs:
            self.sides_ids = kwargs['sides_ids']
        else:
            self.sides_ids = ((),())

        # List of verts, edges, faces all internally indexed
        if 'vert_list' in kwargs and 'face_list' in kwargs:
            if 'CONV' in kwargs:
                if kwargs['CONV'] == True:
                    self.conv_lists(kwargs['vert_list'], kwargs['face_list'])
                elif kwargs['CONV'] == False:
                    self.vert_list = kwargs['vert_list']
                    self.face_list = kwargs['face_list']
                    if 'edge_list' in kwargs:
                        self.edge_list = kwargs['edge_list']
                else:
                    print("Error! CONV flag must be Bool.")
            else:
                print("Error! Missing CONV flag.")
        else:
            self.vert_list = []
            self.face_list = []
            self.edge_list = []

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
    global nm_section_dict
    
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
                sec = NM_section(name=sec_name, sc_id=sec_ids, sc_pts=sec_pts, nghbr_sc_ids=(nghbrs_min,nghbrs_max))
                nm_section_dict[sec_ids] = sec

    return

# Make segements from mesh object
def segment_meshpy(mesh_vert_co_list, mesh_tet_vert_list, sc_tet_list, mesh_tet_nghbr_list, mesh_face_vert_marker_dict, nm_sec, n_seg_plen):
    
    global nm_segment_dict
    
    # Number of tets in this section
    n_tets = len(sc_tet_list)
    
    # Get the distance of projected center points from the endpoint [0] of line seg of all the tets
    tet_proj_dist_list = []
    for tet in sc_tet_list:
        ctr_pt = Vector([0,0,0])
        for vert_id in mesh_tet_vert_list[tet]:
            ctr_pt += Vector(mesh_vert_co_list[vert_id])
        ctr_pt /= 4.0
        tet_proj_dist_list.append((nm_sec.sc_pts[0] - project_pt_line(nm_sec.sc_pts[0],nm_sec.sc_pts[1],ctr_pt)).length)
            
    # Bin them into segments
    # Bin by distance
    min_dist = min(tet_proj_dist_list)
    max_dist = max(tet_proj_dist_list)
    # Number of segments
    n_seg = int((max_dist - min_dist)*n_seg_plen) # int => round down

    # Distance per segment
    d_dist = 1.0/n_seg_plen
    # Make sure there's at least 1 segment
    if n_seg < 1:
        # This section is so short, we can't even put one plane in the middle
        n_seg = 1
    
    # Store the number of segments
    nm_sec.n_seg = n_seg

    # Create a dictionary of segment index to tet indeces
    seg_to_tet_dict = {}
    for i_seg in range(1,n_seg+1):
        seg_to_tet_dict[i_seg] = []
    for i_tet,tet in enumerate(sc_tet_list):
        dist = tet_proj_dist_list[i_tet]
        j = int((dist - min_dist) / d_dist) + 1
        if j>n_seg: # max-distance tet
            j=n_seg
        seg_to_tet_dict[j] += [tet]

    ###
    # Create the segments in the dictionary
    ###

    for i_seg in range(1,n_seg+1):
        seg_name = nm_sec.name + ('_sg_%02d'%i_seg)
        nm_segment_dict[(nm_sec.sc_id[0],nm_sec.sc_id[1],i_seg)] = NM_segment(name=seg_name, sg_id=(nm_sec.sc_id[0],nm_sec.sc_id[1],i_seg))

    ###
    # Find dividing planes
    ###

    ###
    # Segment-segment boundaries
    ###

    # Go through all the segments
    for i_this_seg in list(seg_to_tet_dict.keys())[:-1]:
        # Segement boundary face list for this segment
        sg_brdr_face_list = []
        
        # Tets in this and the next section
        this_tet_ids = seg_to_tet_dict[i_this_seg]
        next_tet_ids = seg_to_tet_dict[i_this_seg+1]
        
        # Go through all the tets in this seg
        for this_tet_id in this_tet_ids:

            # Get this tet's neighbors
            for nghbr_tet_id in mesh_tet_nghbr_list[this_tet_id]:
                # Surfaces ignored
                if nghbr_tet_id != -1:
                    # Is it a boundary tet
                    if nghbr_tet_id in next_tet_ids:
                        # Get the shared verts
                        this_verts = mesh_tet_vert_list[this_tet_id]
                        nghbr_verts = mesh_tet_vert_list[nghbr_tet_id]
                        shared_verts = list(set(this_verts).intersection(set(nghbr_verts)))
                        # Add face
                        sg_brdr_face_list.append(shared_verts)
        
        # Make plane for this segment boundary
        if len(sg_brdr_face_list) > 0:
            sg_brdr_name = nm_sec.name + ('_sg_%02d_B_'%i_this_seg) + nm_sec.name + ('_sg_%02d'%(i_this_seg+1))
            sg_brdr_sides_names = (nm_sec.name + ('_sg_%02d'%i_this_seg),nm_sec.name + ('_sg_%02d'%(i_this_seg+1)))
            sg_brdr_sides_ids = ((nm_sec.sc_id[0],nm_sec.sc_id[1],i_this_seg),(nm_sec.sc_id[0],nm_sec.sc_id[1],i_this_seg+1))
            plane = NM_plane(name=sg_brdr_name, sides_names=sg_brdr_sides_names, sides_ids=sg_brdr_sides_ids, vert_list=mesh_vert_co_list, face_list=sg_brdr_face_list, CONV=True)
            nm_segment_dict[(nm_sec.sc_id[0],nm_sec.sc_id[1],i_this_seg)].planes_sg.append(plane)

    ###
    # Segment's surface boundaries
    ###

    # All face markers possible
    face_markers_possible = nm_sec.face_marker_dict.keys()

    # Go through all the segments
    for i_this_seg in seg_to_tet_dict.keys():
        # Plane surface face list for this segment
        sg_plane_surf_face_list = []
        # Section surface face dict for this segment, for all the different markers
        sg_sc_surf_face_dict = {}
        
        # Tets in this section
        this_tet_ids = seg_to_tet_dict[i_this_seg]

        '''
        # TEMP
        # Make an object out of all the tets
        if nm_sec.name == "sc_02_03":
            tmp_vert_list = []
            tmp_dict = {}
            tmp_face_list = []

            for t in this_tet_ids:
                this_verts = mesh_tet_vert_list[t]
                for v in this_verts:
                    if not v in tmp_dict.keys():
                        tmp_vert_list.append(mesh_vert_co_list[v])
                        tmp_dict[v] = len(tmp_vert_list) - 1
                for ids in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]:
                    tmp = [tmp_dict[this_verts[ids[0]]],tmp_dict[this_verts[ids[1]]],tmp_dict[this_verts[ids[2]]]]
                    tmp.sort()
                    if not tmp in tmp_face_list:
                        tmp_face_list.append(tmp)
            plane = NM_plane(name="TetMesh", sides_names=("a","b"), sides_ids=((1,2,3),(-1)), vert_list=tmp_vert_list, face_list=tmp_face_list, CONV=True)
            plane.make_plane()
        '''

        # Go through all the tets in this seg
        zero_marker_face_list = []
        for this_tet_id in this_tet_ids:

            # Check all the faces
            this_verts = mesh_tet_vert_list[this_tet_id]
            for triplet in list(itertools.combinations(this_verts,3)):
                face_marker = mesh_face_vert_marker_dict[frozenset(triplet)]
                
                # Find the faces that constitute the border of the section
                striplet = sorted(triplet)
                if face_marker in face_markers_possible:
                    # This face is one of the border faces - but what type of face?
                    if face_marker < 0:
                        sg_plane_surf_face_list.append(striplet)
                    elif face_marker in sg_sc_surf_face_dict:
                        sg_sc_surf_face_dict[face_marker].append(striplet)
                    else:
                        sg_sc_surf_face_dict[face_marker] = [striplet]

                # TEMP
                '''
                elif face_marker == 0:
                    # This face is supposedly one of the interior faces
                    # However, the unfortunately the assignments of tets to element attributes by TetGen is not perfect
                    # and does not work as expected. Hence this needs to be checked by counting neighbors

                    # To be an interior face, this face needs to be shared by two tets, i.e. this face needs to occur twice
                    # Keep a running list of all triplets that have occurred to check this
                    if striplet in zero_marker_face_list:
                        del zero_marker_face_list[zero_marker_face_list.index(striplet)]
                    else:
                        zero_marker_face_list.append(striplet)
                '''

        # TEMP
        '''
        # Check that there are no interior faces mal-assigned
        if len(zero_marker_face_list) > 0:
            # These faces are "boundary/surface" faces!
            # Assign them by searching randomly for neighbors
            for zero_f in zero_marker_face_list:
                print("> > Need to assign zero_f = " + str(zero_f))
                # Search through all assigned face markers
                for face_marker, trip_list in sg_sc_surf_face_dict.items(): # NOTE: should NOT need to check sg_plane_surf_face_list
                    # Flatten the list
                    trip_flat = [item for sublist in trip_list for item in sublist]
                    if zero_f[0] in trip_flat and zero_f[1] in trip_flat:
                        # Just take this surface face marker
                        sg_sc_surf_face_dict[face_marker].append(zero_f)
                        print("Assigned to face marker: " + str(face_marker))
                        break # Stop!
        '''

        # TEMP
        '''
        # Make an object out of all the faces that constitute the surface
        if nm_sec.name == "sc_02_03":
            tmp_face_list = sg_plane_surf_face_list
            for l in sg_sc_surf_face_dict.values():
                tmp_face_list += l
            # Vert list
            tmp_dict = {}
            tmp_vert_list = []
            for f in tmp_face_list:
                for v in f:
                    if not v in tmp_dict.keys():
                        tmp_vert_list.append(mesh_vert_co_list[v])
                        tmp_dict[v] = len(tmp_vert_list) - 1
            tmp_face_list_2 = []
            for f in tmp_face_list:
                tmp_face_list_2.append([tmp_dict[f[0]],tmp_dict[f[1]],tmp_dict[f[2]]])
            plane = NM_plane(name="Border", sides_names=("a","b"), sides_ids=((1,2,3),(-1)), vert_list=tmp_vert_list, face_list=tmp_face_list_2, CONV=True)
            plane.make_plane()
        '''

        # Make the surface plane for this segment
        if len(sg_plane_surf_face_list) > 0:
            sg_surf_name = nm_sec.name + ('_sg_%02d_B_surf'%i_this_seg)
            sg_surf_sides_names = (nm_sec.name + ('_sg_%02d'%i_this_seg),'surf')
            sg_surf_sides_ids = ((nm_sec.sc_id[0],nm_sec.sc_id[1],i_this_seg),(-1))
            plane = NM_plane(name=sg_surf_name, sides_names=sg_surf_sides_names, sides_ids=sg_surf_sides_ids, vert_list=mesh_vert_co_list, face_list=sg_plane_surf_face_list, CONV=True)
            nm_segment_dict[(nm_sec.sc_id[0],nm_sec.sc_id[1],i_this_seg)].plane_surf = plane

        # Make the section-surface plane for this segment
        if len(sg_sc_surf_face_dict.keys()) > 0:
            for face_marker in sg_sc_surf_face_dict.keys():
                sg_sc_surf_face_list = sg_sc_surf_face_dict[face_marker]
                # The hard part: wtf am i bordering?
                brdring_plane_idxs = nm_sec.face_marker_dict[face_marker]
                brdring_plane = nm_sec.planes_sc_brdrs[brdring_plane_idxs[0]][brdring_plane_idxs[1]]
                sides_ids = brdring_plane.sides_ids
                if sides_ids[0] == nm_sec.sc_id:
                    brdring_id = sides_ids[1]
                elif sides_ids[1] == nm_sec.sc_id:
                    brdring_id = sides_ids[0]
                else:
                    print("Error! Something went wrong with your planes.")
                
                min_id = min(brdring_id,nm_sec.sc_id)
                max_id = max(brdring_id,nm_sec.sc_id)
                if min_id == nm_sec.sc_id:
                    sg_surf_name = 'sc_%02d_%02d_sg_%02d_B_'%(min_id[0],min_id[1],i_this_seg) + 'sc_%02d_%02d'%(max_id[0],max_id[1])
                    sg_surf_sides_names = ('sc_%02d_%02d_sg_%02d'%(min_id[0],min_id[1],i_this_seg),'sc_%02d_%02d'%(max_id[0],max_id[1]))
                    sg_surf_sides_ids = ((min_id[0],min_id[1],i_this_seg),(max_id[0],max_id[1]))
                else:
                    sg_surf_name = 'sc_%02d_%02d_B_'%(min_id[0],min_id[1]) + 'sc_%02d_%02d_sg_%02d'%(max_id[0],max_id[1],i_this_seg)
                    sg_surf_sides_names = ('sc_%02d_%02d'%(min_id[0],min_id[1]),'sc_%02d_%02d_sg_%02d'%(max_id[0],max_id[1],i_this_seg))
                    sg_surf_sides_ids = ((min_id[0],min_id[1]),(max_id[0],max_id[1],i_this_seg))

                plane = NM_plane(name=sg_surf_name, sides_names=sg_surf_sides_names, sides_ids=sg_surf_sides_ids, vert_list=mesh_vert_co_list, face_list=sg_sc_surf_face_list, CONV=True)
                nm_segment_dict[(nm_sec.sc_id[0],nm_sec.sc_id[1],i_this_seg)].planes_sc.append(plane)

                '''
                # TEMP
                if min_id[0] == 2 and min_id[1] == 3 and max_id[0] == 2 and max_id[1] == 6:
                    print(">>>> Making: " + str(sg_surf_name))
                    plane.make_plane()
                '''

    # Fin!
    return

# New object from MeshPy
def new_obj_meshpy(obj_name, mesh_built):
    
    mesh_new = bpy.data.meshes.new(obj_name+"_mesh")
    
    vert_list = list(mesh_built.points)
        
    edge_list = list(mesh_built.edges)
    
    face_list = list(mesh_built.faces)
    
    # Build the object
    mesh_new.from_pydata(vert_list,edge_list,face_list)
    
    # Update
    mesh_new.validate(verbose=False) # Important! and i dont know why
    mesh_new.update()

    # Object
    obj_new = bpy.data.objects.new(obj_name,mesh_new)

    # Something
    scene = bpy.context.scene
    scene.objects.link(obj_new)

# Main

def f_compartmentize(context, swc_filepath):

    print("> Running: f_compartmentize")

    # Get the working dir from the filepath
    global DIR
    swc_filename = os.path.basename(swc_filepath)
    DIR = swc_filepath[:-len(swc_filename)]
    print("Working directory: " + str(DIR))

    # Time
    t_st = []
    t_st.append(time.time())

    global nm_section_dict
    global nm_segment_dict
    
    nm_section_dict = {}
    nm_segment_dict = {}
    
    # Get data from the SWC file
    get_connections(swc_filepath)

    # Get the active object
    ob_list = context.selected_objects

    if len(ob_list) != 1:
        return
    else:
        ob = ob_list[0]
    
    # All vertices
    ob_vert_list = [tuple(item.co) for item in ob.data.vertices]

    # Make sure everything is de-selected before we start
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')

    # Get the regions of the cube object
    reg_list = ob.mcell.regions.region_list

    # Go over all sections, get the border edges
    sc_edge_dict = {}
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
            
            # Get the faces
            bpy.ops.object.mode_set(mode='EDIT')
            sec.select_region_faces(context)
            
            # Store the surface as a new plane
            bpy.ops.object.mode_set(mode='OBJECT')
            surf_face_list = [i.vertices for i in ob.data.polygons if i.select]
            surf_name = 'sc_%02d_%02d'%(pt_min,pt_max)
            surf_sides_names = (surf_name,-1)
            surf_sides_ids = ((pt_min,pt_max),(-1))
            # Make the plane
            plane_surf = NM_plane(name=surf_name, sides_names=surf_sides_names, sides_ids=surf_sides_ids, vert_list=ob_vert_list, face_list=surf_face_list, CONV=True)
            # Store the plane as the surface for this section
            nm_section_dict[(pt_min,pt_max)].plane_surf = plane_surf
            
            # Get the edges
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.region_to_loop()
            bpy.ops.object.mode_set(mode='OBJECT')
            sc_edge_dict[(pt_min,pt_max)] = [i.index for i in ob.data.edges if i.select]

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

    # Go through all of the sections
    for this_sc_id in nm_section_dict.keys():
        this_sc = nm_section_dict[this_sc_id]
        
        # Border
        this_sc_edge = sc_edge_dict[this_sc_id]
        
        # Go through both of end verts
        for i_end,nghbr_sc_ids in enumerate(this_sc.nghbr_sc_ids):
            
            # Go through all of the neighboring sections
            for nghbr_sc_id in nghbr_sc_ids:
                
                # Don't double count
                if nghbr_sc_id > this_sc_id:
                    
                    # Border
                    nghbr_sc_edge = sc_edge_dict[nghbr_sc_id]

                    # Get the elements they share
                    edge_share = list(set(this_sc_edge).intersection(set(nghbr_sc_edge)))

                    if len(edge_share) > 0: # Do they share any?
                        
                        # Vertices
                        plane_vert_list = ob_vert_list + [tuple(this_sc.sc_pts[i_end])]
                        
                        # Index of the ctr pt
                        ctr_pt_id = len(plane_vert_list)-1
                        
                        # Convert the edges to vertices
                        verts_share = [[ob.data.edges[edge].vertices[0],ob.data.edges[edge].vertices[1]] for edge in edge_share]

                        # Faces
                        plane_face_list = [[ctr_pt_id, item[0], item[1]] for item in verts_share]

                        # Get the actual neighboring section
                        nghbr_sc = nm_section_dict[nghbr_sc_id]

                        # Make a new plane
                        name_plane = this_sc.name + "_B_" + nghbr_sc.name
                        plane = NM_plane(name=name_plane,sides_names=(this_sc.name,nghbr_sc.name),sides_ids=(this_sc_id,nghbr_sc_id),vert_list=plane_vert_list,face_list=plane_face_list,CONV=True)
                        
                        # Add the plane to both section boundaries
                        this_sc.planes_sc_brdrs[i_end].append(plane)
                        if nghbr_sc_id[0] == this_sc_id[i_end]:
                            nghbr_sc.planes_sc_brdrs[0].append(plane)
                        elif nghbr_sc_id[1] == this_sc_id[i_end]:
                            nghbr_sc.planes_sc_brdrs[1].append(plane)
                        else:
                            print("Something went wrong!")

    # # # Free memory
    del ob

    # Time
    t_st.append(time.time())
    print("> Made section boundaries: time: " + str(t_st[-1]-t_st[-2]))

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

    for sec_id in list(nm_section_dict.keys()):
        
        # if sec_id in sc_tet_dict: # Why should this not be the case? This occurs for some reason....
        
        sec = nm_section_dict[sec_id]

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
                                        plane = NM_plane(name=p_name, sides_names=p_sides_names, sides_ids=p_sides_ids, vert_list=vert_list, face_list=face_list, CONV=True)
                                        
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
    nm_surface_plane = NM_plane(name="NM_Surface",vert_list = plane_surf_vert_list, face_list = plane_surf_face_list, CONV=True)

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

    nm_seg_plane = NM_plane(name="NM_Segment",vert_list = tmp_vert_list, face_list = tmp_face_list, CONV=True)
    '''

    # Make the plane
    nm_seg_plane = NM_plane(name="NM_Segment",vert_list = seg_surf_vert_list, face_list = seg_surf_face_list, CONV=True)

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
            sc_pts = nm_section_dict[(id_min_min,id_min_max)].sc_pts
            vec_check = sc_pts[1]-sc_pts[0]
        else:
            # Section boundary
            all_ids = [id_min_min, id_min_max, id_max_min, id_max_max]
            shared_id = list(set([item for item in all_ids if all_ids.count(item) == 2]))[0]
            uniq_ids = [item for item in all_ids if all_ids.count(item) == 1]
            max_uniq_id = max(uniq_ids)
            min_id = min(shared_id,max_uniq_id)
            max_id = max(shared_id,max_uniq_id)
            sc_pts = nm_section_dict[(min_id,max_id)].sc_pts
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








