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

# Function to construct a neighbour list from a list of face vertices
def construct_neighbor_list(face_list):
    neigh_list = []
    for f in face_list:
        neigh_list.append([])

    print("Constructing neighbor list....")

    # Go through each face
    for i_f, f in enumerate(face_list):
        
        if i_f % 100 == 0:
            print(str(i_f) + " / " + str(len(face_list)))
        
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
                    neigh_list[i_f].append(i_fo+i_f+1)
                    neigh_list[i_fo+i_f+1].append(i_f)
                    break

    return neigh_list

# Main

def f_check_bordering_surfaces(context, swc_filepath):

    print("> Running: f_check_bordering_surfaces")

    ###
    # Read the connections based on SWC file
    ###
    
    # Get data from the SWC file
    pt_connect = get_connections(swc_filepath)

    ###
    # Check connections in the file
    ###
    
    # Get the object
    ob_list = context.selected_objects

    if len(ob_list) == 1:
        ob = ob_list[0]
        
        # List of faces defined by the vertices
        face_list = [item.vertices for item in ob.data.polygons]
        
        # Construct dictionary of neighbors of all faces
        neigh_list = construct_neighbor_list(face_list)
            
        # Get the mcell regions
        reg_list = ob.mcell.regions.region_list
        
        # Construct dictionary based on names
        reg_dict = {}
        for reg in reg_list:
            # Id of this region
            reg_id = (int(reg.name[3:5]),int(reg.name[6:8]))
            # Get region's faces
            reg_face_list = list(reg.get_region_faces(ob.data))
            # Store
            reg_dict[reg_id] = reg_face_list

        print("Checking MCell region's face assignments")

        # Store problematic faces
        f_problem_list = []

        # Go through all mcell regions
        for reg_id, reg_face_list in reg_dict.items():
            
            print("Checking region: " + str(reg_id))
            
            # Construct all the neighboring face lists allowed by the SWC file
            neighs_allowed = []
            id_pairs_counted = []
            for i_conn in [0,1]:
                for o_conn in pt_connect[reg_id[i_conn]-1]:
                    id_pair = (min(reg_id[i_conn],o_conn),max(reg_id[i_conn],o_conn))
                    if not id_pair in id_pairs_counted:
                        neighs_allowed += reg_dict[id_pair]
                        id_pairs_counted.append(id_pair)

            # Go through every face
            for f in reg_face_list:
            
                # Check neighbors
                neighbors = neigh_list[f]
                
                # Check if the neighbors are allowed based on the SWC file
                for neigh in neighbors:
                    # Not in the same region or a neighboring region
                    if not (neigh in reg_face_list or neigh in neighs_allowed):
                        print("---   ---   ---")
                        print("Error!")
                        print("Problematic face #: " + str(f) + " is in region " + str(reg_id))
                        # Find what section it's in instead
                        for reg_o_id, reg_o_face_list in reg_dict.items():
                            if neigh in reg_o_face_list:
                                print("Borders face #: " + str(neigh) + " is in region " + str(reg_o_id))
                                print("This is not allowed by the SWC file!")
                                print("Region: " + str(reg_id) + " can only border: ")
                                for i_conn in [0,1]:
                                    for o_conn in pt_connect[reg_id[i_conn]-1]:
                                        id_pair = (min(reg_id[i_conn],o_conn),max(reg_id[i_conn],o_conn))
                                        if id_pair != reg_id:
                                            print(id_pair)
                                print("---   ---   ---")

                                # Store the problematic face ids
                                f_problem_list.append(f)
                                f_problem_list.append(neigh)

                                # Stop for this face (but keep looking for other faces in violation)
                                break

        # Select problematic faces if they exist
        if len(f_problem_list) > 0:

            # First set the mode
            bpy.ops.object.mode_set(mode='EDIT')
            # Face selection mode
            context.tool_settings.mesh_select_mode = (False, False, True)
            # Deselect all
            bpy.ops.mesh.select_all(action='DESELECT')
            # Must be in object mode to select faces!
            bpy.ops.object.mode_set(mode='OBJECT')

            # Select the problematic faces
            for f in f_problem_list:
                ob.data.polygons[f].select = True

            # Return in edit mode
            bpy.ops.object.mode_set(mode='EDIT')

    print("> Finished: f_check_bordering_surfaces")


