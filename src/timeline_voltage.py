import bpy

import math

import struct

import os
import sys

import cellblender

# Function to read a single data file
def read_voltage_data(fname, v_dict):
	f = open(fname, 'r')

	for i, line in enumerate(f):
		line_s = line.split();
		if len(line_s) == 4:
			key = (int(line_s[0]),int(line_s[1]),int(line_s[2]))
			v_dict[key] = float(line_s[3])

	f.close()

	return


'''
# Set colors of the molecule glyphs
def set_mol_cols():
	
	# Colors
	col_dict = {}
	col_off = (0.9,0.9,0.9)
	col_dict['s00'] = col_off
	col_dict['s10'] = col_off
	col_dict['s20'] = col_off
	col_dict['s30'] = col_off
	col_dict['s01'] = col_off
	col_dict['s11'] = col_off
	col_dict['s21'] = col_off
	col_dict['s31'] = (0.0,1.0,0.0)

	# Set colors
	mol_par_obj = bpy.data.objects.get("molecules")
	
	# Check that it exists
	if mol_par_obj:
		# Get children
		for mol_obj in mol_par_obj.children:
			shape_obj = mol_obj.children[0]
			
			# Are there any materials?
			if len(shape_obj.data.materials)>0:
				
				# Get the material
				mat = shape_obj.data.materials[0]
				mol_name = mn[4:7]
				mat.diffuse_color = col_dict[mol_name]
'''

# Main

def f_timeline_voltage(scene, frame):

	# Print to console to delimit
	# print("> Running f_timeline_voltage")

	# Go through all the objects
	for mesh_obj in scene.nrnlauncher.mesh_obj_list:

		# Check that it has a voltage list
		v_dir = mesh_obj.v_dir
		v_zero_pad = mesh_obj.v_zero_pad
		v_n_files = mesh_obj.v_n_files
		if v_dir != "":

			# The voltage file no we want to read
			v_file_no = int(frame) % v_n_files

			# Get the object
			ob = bpy.data.objects[mesh_obj.name]

			# Read the data for this object at this frame
			v_dict = {}
			read_voltage_data((v_dir+"v_%0"+str(v_zero_pad)+"d.txt")%v_file_no, v_dict)

			# Check that there is some voltage data to proceed
			if not len(v_dict) > 0:
				raise SystemError("No voltage data read!")

			# Max and min voltages for colors
			min_v = -80
			max_v = 20

			# Go through all materials        
			# mat_list = [item.material for item in bpy.context.object.material_slots]
			for mat in ob.data.materials:
				mn = mat.name
				if len(mn) == 23 and mn[-9:] == "_material":
					reg_id = (int(mn[3:5]),int(mn[6:8]),int(mn[12:14]))
					if reg_id in v_dict:
						# Change color
						frac = min(1.0,max(0.0,(v_dict[reg_id] - min_v)/abs(max_v - min_v)))
						mat.diffuse_color = (frac,0.0,1.0-frac)

	# print("> Finished f_timeline_voltage")

