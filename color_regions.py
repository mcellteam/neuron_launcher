import bpy

import math

import struct

import os
import sys

import cellblender

import random

# Print to console to delimit
print("--- Running timeline_voltage.py ---")

# Take regions in MCell and make materials for each
def f_mcell_reg_to_mat(context):
    
    print("> Running: f_mcell_reg_to_mat")

    # Get the active object
    ob_list = context.selected_objects

    if len(ob_list) == 1:
        ob = ob_list[0]

        # Get the regions of the cube object
        reg_list = ob.mcell.regions.region_list

        # Clear all materials on the object
        bpy.ops.object.mode_set(mode='OBJECT')
        for i in range(0,context.object.material_slots.__len__()):
            context.object.active_material_index = 1
            bpy.ops.object.material_slot_remove()
        ob.data.materials.clear()

        # Assign a material to each region
        for i_reg in reg_list:
            
            if len(i_reg.name) == 8 and i_reg.name[0:2] == 'sc':
                
                # Deselect all objects currently selected
                # bpy.ops.object.select_all(action='DESELECT')
                
                # Select all faces in this region > Edit mode
                i_reg.select_region_faces(bpy.context)

                # Add a material slot
                bpy.ops.object.material_slot_add()

                # Assign a material to the last slot
                context.object.material_slots[bpy.context.object.material_slots.__len__() - 1].material = bpy.data.materials.new(name=i_reg.name+"_material")
                
                # The name automatically gets something appended to it... change it back to the name we want
                context.object.material_slots[context.object.material_slots.__len__() - 1].material.name = i_reg.name+"_material"

                #Assign the material on the selected vertices
                bpy.ops.object.material_slot_assign()
                
                # Deselect all vertices
                bpy.ops.mesh.select_all(action='DESELECT')

    print("> Finished: f_mcell_reg_to_mat")

# Color regions blue
def f_color_regions(context):

    print("> Running: f_color_regions")

    # Get the active object
    ob_list = context.selected_objects

    if len(ob_list) == 1:
        ob = ob_list[0]

        # Go through all materials
        mat_list = [item.material for item in context.object.material_slots]
        for mat in mat_list:
            i = random.random()
            mat.diffuse_color = (0.0,i,1.0-i)
            mat.use_transparency = True
            mat.transparency_method = 'MASK'
            mat.alpha = 0.4

    print("> Finished: f_color_regions")





