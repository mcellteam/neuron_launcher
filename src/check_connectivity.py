import bpy, bmesh
from mathutils import Vector
import math

import sys

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

    return pt_connect

# Main

def f_check_connectivity(context, swc_filepath):

    print("> Running: f_check_connectivity")

    ###
    # Read the connections based on SWC file
    ###
    
    # Get data from the SWC file
    pt_connect = get_connections(swc_filepath)

    # Get the object
    ob_list = context.selected_objects

    if len(ob_list) != 1:
        raise TypeError("Please select only one object.")

    ob = ob_list[0]

    # Get the mcell regions
    reg_list = ob.mcell.regions.region_list
    
    # Dictionary of bordering region idxs
    reg_brdr_dict = {}
    for reg_id, reg in enumerate(reg_list):
        reg_brdr_dict[reg_id] = []

    # Go through every region
    for reg_id, reg in enumerate(reg_list):

        print("Checking region " + str(reg_id) + " / " + str(len(reg_list)-1) + " name: " + str(reg.name) + " ....")

        # Deselect all
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT')

        # Select the region's faces
        reg.select_region_faces(context)

        # Convert to edge loop
        bpy.ops.mesh.region_to_loop()

        # Get the edges
        bpy.ops.object.mode_set(mode='OBJECT')
        edge_loop = [i.index for i in ob.data.edges if i.select]

        # Convert to vertices
        verts_loop = [tuple(ob.data.edges[edge].vertices) for edge in edge_loop]

        # Go through all faces in the object, find the faces that share these edges
        f_list = []
        for f in ob.data.polygons:
            es = list(f.edge_keys)
            for e in es:
                if e in verts_loop:
                    f_list.append(f.index)
                    break

            # Maximum size (found all the faces)
            if len(f_list) == 2*len(edge_loop):
                break

        # Get the MCell region of each face
        reg_brdr_list = []
        for o_reg_id, o_reg in enumerate(reg_list[reg_id+1:]):
            o_f_list = o_reg.get_region_faces(ob.data)
            for f in f_list:
                if f in o_f_list:
                    # This region borders
                    reg_brdr_dict[reg_id].append(reg_id+1+o_reg_id)
                    reg_brdr_dict[reg_id+1+o_reg_id].append(reg_id)
                    break

        # Convert the reg ids into names
        reg_brdr_names = []
        for o_reg_id in reg_brdr_dict[reg_id]:
            o_reg = reg_list[o_reg_id]
            if not o_reg.name in reg_brdr_names:
                reg_brdr_names.append(o_reg.name)

        # Construct list of allowed region names based on SWC file connectivity
        pt1 = int(reg.name[3:5])
        pt2 = int(reg.name[6:8])
        reg_allowed_names = []
        for i_pt in [pt1,pt2]:
            for j_pt in pt_connect[i_pt-1]:
                allowed_name = "sc_%02d_%02d" % (min(i_pt,j_pt),max(i_pt,j_pt))
                if not allowed_name in reg_allowed_names and not reg.name == allowed_name:
                    reg_allowed_names.append(allowed_name)

        # Go through the brdr regions and check if they're allowed
        bad_regs = []
        for brdr_name in reg_brdr_names:
            if not brdr_name in reg_allowed_names:
                bad_regs.append(brdr_name)

        # Go through the brdr regions and check if one is missed
        missed_regs = []
        for allowed_name in reg_allowed_names:
            if not allowed_name in reg_brdr_names:
                missed_regs.append(allowed_name)

        # Print and inform!
        if len(bad_regs) > 0:
            print("")
            print("ERROR: Region: " + str(reg.name))
            print("illegally borders the following:")
            for brdr_name in bad_regs:
                print("-N- " + str(brdr_name))
            print("only the following are allowed:")
            for allowed_name in reg_allowed_names:
                print("-Y- " + str(allowed_name))
            print("")
        if len(missed_regs) > 0:
            print("")
            print("WARNING: Region: " + str(reg.name))
            print("fails to border the following (this may or may not be ok):")
            for missed_name in missed_regs:
                print("-F- " + str(missed_name))
            print("")
        if len(bad_regs) == 0 and len(missed_regs) == 0:
            print("- Region borders are OK!")

    # Return in edit mode
    bpy.ops.object.mode_set(mode='EDIT')
    context.tool_settings.mesh_select_mode = (False, False, True)

    print("> Finished: f_check_connectivity")



