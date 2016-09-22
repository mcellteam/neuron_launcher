import bpy, bmesh
from mathutils import Vector
import math

import collections

import numpy as np

import functools

import itertools

import os

# Time
import time

# Bisect
import bisect

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

    swc_vert_list = [Vector(item[2:5]) for item in swc_data]

    return pt_connect, swc_vert_list

# Main

def f_explode(context, swc_filepath, exp_scale_factor):

    print("> Running: f_explode")

    # Time
    t_st = []
    t_st.append(time.time())

    # Get data from the SWC file
    pt_connect, swc_vert_list = get_connections(swc_filepath)

    #######################################################
    #######################################################
    # Go through each section, figure out the displacement vector
    #######################################################
    #######################################################

    # First, get a list of all the center points of each section
    ctr_dict = {}
    for i_pt, conn_pts in enumerate(pt_connect):
        for i_conn_pt_p in conn_pts:
            i_conn_pt = i_conn_pt_p - 1
            if i_conn_pt > i_pt: # Don't double count sections
                co_1 = swc_vert_list[i_pt]
                co_2 = swc_vert_list[i_conn_pt]
                # Center
                ctr = co_1 + 0.5*(co_2 - co_1)
                ctr_dict[(i_pt+1,i_conn_pt+1)] = ctr

    # Compute the global center
    ctr_global = Vector([0,0,0])
    for c in ctr_dict.values():
        ctr_global += c
    ctr_global /= len(ctr_dict)
    print(ctr_global)

    # Subtract off the global center
    for key, val in ctr_dict.items():
        ctr_dict[key] = val - ctr_global

    # Scale all the points and store the displacement vectors
    disp_dict = {}
    for key,val in ctr_dict.items():
        disp_dict[key] = exp_scale_factor * val - val

    # Go through all objects
    for ob0 in context.scene.objects:
        oname = ob0.name

        # Check the name that it is a section we wish to displace
        if len(oname) == 20 and oname[0:3] == "sc_" and oname[-11:] == "compartment":

            print("Exploding: " + str(oname))

            id1 = int(oname[3:5])
            id2 = int(oname[6:8])

            disp = disp_dict[(id1,id2)]

            # Get the Bmesh
            if ob0.data.is_editmode:
                bm = bmesh.from_edit_mesh(ob0.data)
            else:
                bm = bmesh.new()
                bm.from_mesh(ob0.data)
            
            # Translate
            bmesh.ops.translate(bm, verts=bm.verts, vec = tuple(disp))

            # Update object
            if bm.is_wrapped:
                bmesh.update_edit_mesh(ob0.data)
            else:
                bm.to_mesh(ob0.data)
                ob0.data.update()

    print("> Finished: f_explode")








