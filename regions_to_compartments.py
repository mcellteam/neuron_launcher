import bpy, bmesh
from mathutils import Vector
import math

import collections

import functools

import itertools

import os

# To add objects to MCell
from cellblender.cellblender_utils import preserve_selection_use_operator

# Main

def f_regions_to_compartments(context):

	print("> Running: f_regions_to_compartments")

	############################
	############################
	# Identify which object is the surface and which is the segment
	############################
	############################

	# Get the active objects
	ob_list = context.selected_objects

	if len(ob_list) != 2:
		raise SystemError("Please select two objects: the surface object (with surface patches) and the segment object (with dividing planes)")

	# Determine which object is the surface and which the segment
	reg_list = ob_list[0].mcell.regions.region_list

	if len(reg_list) == 0:
		raise SystemError("One or both objects have no MCell regions defined.")

	sample_name = reg_list[0].name

	if len(sample_name) >= 6 and sample_name[-6:] == "B_surf":
		# This is the surface object
		ob_surf = ob_list[0]
		ob_seg = ob_list[1]
	else:
		# This is the segment object
		ob_surf = ob_list[1]
		ob_seg = ob_list[0]

	reg_list_surf = ob_surf.mcell.regions.region_list
	reg_list_seg = ob_seg.mcell.regions.region_list

	# Check that the lengths of the names in each are all the same
	len_surf = len(reg_list_surf[0].name)
	if len(reg_list_surf) > 1:
		for reg in reg_list_surf[1:]:
			if len(reg.name) != len_surf:
				raise SystemError("The surface object has regions of variable name lengths, i.e. a naming error.")
	len_seg = len(reg_list_seg[0].name)
	if len(reg_list_seg) > 1:
		for reg in reg_list_seg[1:]:
			if len(reg.name) != len_seg:
				raise SystemError("The segment object has regions of variable name lengths, i.e. a naming error.")

	# Check that the lengths of the names are compatible. Either:
	# - Sections only; names are:
	# - - Surface: sc_##_##_B_surf = length 15
	# - - Segment: sc_##_##_B_sc_##_## = length 19
	# Segments: names are
	# - - Surface: sc_##_##_sg_##_B_surf = length 21
	# - - Segment: sc_##_##_sg_##_B_sc_##_##_sg_## = length 31
	if not ((len_surf == 15 and len_seg == 19) or (len_surf == 21 and len_seg == 31)):
		raise SystemError("The surface and segment objects have incompatible region name lengths.")

	############################
	############################
	# Make the compartment objects
	############################
	############################

	# Turn on edit mode
	bpy.ops.object.mode_set(mode='EDIT')

	# Get a dictionary of all materials in the surface object
	context.scene.objects.active = ob_surf
	mat_list = [item.material for item in context.object.material_slots]
	mat_dict = {}
	for mat in mat_list:
		mat_name = mat.name
		if len_seg == 19:
			sc_id = (int(mat_name[3:5]),int(mat_name[6:8]))
		elif len_seg == 31:
			sc_id = (int(mat_name[3:5]),int(mat_name[6:8]),int(mat_name[12:14]))
		mat_dict[sc_id] = mat

	# Go through all the surface regions
	for reg_surf in reg_list_surf:
		name_surf = reg_surf.name
		if len_seg == 19:
			sc_id = (int(name_surf[3:5]),int(name_surf[6:8]))
		elif len_seg == 31:
			sc_id = (int(name_surf[3:5]),int(name_surf[6:8]),int(name_surf[12:14]))

		print("Making object for compartment: " + str(sc_id))

		# Reset the vertex, edge, face lists
		v_new_list = []
		v_id_dict = {} # From the surface idxs to the new vert list
		e_new_list = []
		f_new_list = []

		# Get the faces in this region
		f_reg_surf_list = reg_surf.get_region_faces(ob_surf.data)

		# Add to the lists

		# Vertices
		f_v_list = [ob_surf.data.polygons[item].vertices for item in f_reg_surf_list]
		v_unique_list = [item for sublist in f_v_list for item in sublist]
		v_unique_list = list(set(v_unique_list))
		for v in v_unique_list:
			v_new_list.append(tuple(ob_surf.data.vertices[v].co))
			v_id_dict[v] = len(v_new_list) - 1
		# Edges, faces
		for fs in f_v_list:
			f_new_list.append([v_id_dict[v] for v in fs])
			e_new_list.append([v_id_dict[fs[0]], v_id_dict[fs[1]]])
			e_new_list.append([v_id_dict[fs[0]], v_id_dict[fs[2]]])
			e_new_list.append([v_id_dict[fs[1]], v_id_dict[fs[2]]])

		# Reset the id dict; now it's the dictionary from the segment idxs to the new vert list
		v_id_dict = {}

		# Go through all segments
		for reg_seg in reg_list_seg:
			if len_seg == 19:
				sc_id_1 = (int(reg_seg.name[3:5]),int(reg_seg.name[6:8]))
				sc_id_2 = (int(reg_seg.name[14:16]),int(reg_seg.name[17:19]))
			elif len_seg == 31:
				sc_id_1 = (int(reg_seg.name[3:5]),int(reg_seg.name[6:8]),int(reg_seg.name[12:14]))
				sc_id_2 = (int(reg_seg.name[20:22]),int(reg_seg.name[23:25]),int(reg_seg.name[29:31]))

			# Is this a boundary we care about?
			if sc_id_1 == sc_id or sc_id_2 == sc_id:
				# Yes it is...

				# Get the faces in this region
				f_reg_seg_list = reg_seg.get_region_faces(ob_seg.data)

				# Add to the lists

				# Vertices
				f_v_list = [ob_seg.data.polygons[item].vertices for item in f_reg_seg_list]
				v_unique_list = [item for sublist in f_v_list for item in sublist]
				v_unique_list = list(set(v_unique_list))
				v_unique_co_list = [tuple(ob_seg.data.vertices[v].co) for v in v_unique_list]
				# Check if it already exists
				for i_v, vco in enumerate(v_unique_co_list):
					v = v_unique_list[i_v]
					if not v in v_id_dict:
						# Search if the vertex already exists from the surface
						try:
							idx = v_new_list.index(vco)
							v_id_dict[v] = idx 
						# Append it as a new vertex
						except ValueError:
							v_new_list.append(vco)
							v_id_dict[v] = len(v_new_list) - 1

				# Edges, faces
				for fs in f_v_list:
					f_new_list.append([v_id_dict[v] for v in fs])
					e_new_list.append([v_id_dict[fs[0]], v_id_dict[fs[1]]])
					e_new_list.append([v_id_dict[fs[0]], v_id_dict[fs[2]]])
					e_new_list.append([v_id_dict[fs[1]], v_id_dict[fs[2]]])

		###################################
		###################################
		# Finally: make the object
		###################################
		###################################

		if len_seg == 19:
			obj_new_name =  "sc_%02d_%02d_compartment" % sc_id
		elif len_seg == 31:
			obj_new_name = "sc_%02d_%02d_sg_%02d_compartment" % sc_id
		mesh_new = bpy.data.meshes.new(obj_new_name + "_mesh")
		mesh_new.from_pydata(v_new_list,e_new_list,f_new_list)
		# Validate and update
		mesh_new.validate(verbose=False) # Important! and i dont know why
		mesh_new.update()
		# Overwrite existing object
		obj_old = bpy.data.objects.get(obj_new_name)
		if obj_old:
			context.scene.objects.unlink(obj_old)
			bpy.data.objects.remove(obj_old)
		# New one
		obj_new = bpy.data.objects.new(obj_new_name,mesh_new)
		# Link
		context.scene.objects.link(obj_new)

		###################################
		###################################
		# Add the material to the object
		###################################
		###################################

		# Make sure the new object is the active
		context.scene.objects.active = obj_new

		# Select all vertices
		for f in obj_new.data.polygons:
			f.select = True

		# Add a material slot
		bpy.ops.object.material_slot_add()

		# Assign a material to the last slot
		context.object.material_slots[bpy.context.object.material_slots.__len__() - 1].material = mat_dict[sc_id]

		# Assign the material on the selected vertices
		bpy.ops.object.material_slot_assign()
		

	print("> Finished: f_regions_to_compartments")







