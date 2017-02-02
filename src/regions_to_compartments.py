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

	# Find a name of a compartment to test the length
	for reg in reg_list_surf:
		if ( len(reg.name) == 15 and reg.name[0:3] == "sc_" ) or ( len(reg.name) == 21 and reg.name[0:3] == "sc_" ):
			len_surf = len(reg.name)
			break
	for reg in reg_list_seg:
		if ( len(reg.name) == 19 and reg.name[0:3] == "sc_" ) or ( len(reg.name) == 31 and reg.name[0:3] == "sc_" ):
			len_seg = len(reg.name)
			break

	# Check that the lengths of the names are compatible. Either:
	# - Sections only; names are:
	# - - Surface: sc_##_##_B_surf = length 15
	# - - Segment: sc_##_##_B_sc_##_## = length 19
	# Segments: names are
	# - - Surface: sc_##_##_sg_##_B_surf = length 21
	# - - Segment: sc_##_##_sg_##_B_sc_##_##_sg_## = length 31
	if not ((len_surf == 15 and len_seg == 19) or (len_surf == 21 and len_seg == 31)):
		raise SystemError("The surface and segment objects have incompatible compartment region name lengths.")

	############################
	############################
	# Go through, find the face lists of all materials on both objects
	############################
	############################

	print("Making lists of all materials faces....")

	# Material name to face indexes
	mat_surf_dict = {}
	mat_seg_dict = {}

	# Deselect both objects
	ob_surf.select = False
	ob_seg.select = False

	# Select the surface object
	ob_surf.select = True
	context.scene.objects.active = ob_surf

	# Start with the first index
	ob_surf.active_material_index = 0

	# Go through all materials
	while ob_surf.active_material_index < len(context.object.material_slots):
		# Select
		bpy.ops.object.mode_set(mode='EDIT')
		bpy.ops.mesh.select_all(action='DESELECT')
		bpy.ops.object.material_slot_select()
		bpy.ops.object.mode_set(mode='OBJECT')

		# Count
		f_sel = []
		for f in ob_surf.data.polygons:
			if f.select == True:
				f_sel.append(f.index)

		# Store
		mat = context.object.material_slots[ob_surf.active_material_index].material
		mat_surf_dict[mat.name] = (f_sel,mat)

		# Next!
		ob_surf.active_material_index += 1

	# Deselect this object
	ob_surf.select = False

	# Select the segment object
	ob_seg.select = True
	context.scene.objects.active = ob_seg

	# Start with the first index
	ob_seg.active_material_index = 0

	# Go through all materials
	while ob_seg.active_material_index < len(context.object.material_slots):
		# Select
		bpy.ops.object.mode_set(mode='EDIT')
		bpy.ops.mesh.select_all(action='DESELECT')
		bpy.ops.object.material_slot_select()
		bpy.ops.object.mode_set(mode='OBJECT')

		# Count
		f_sel = []
		for f in ob_seg.data.polygons:
			if f.select == True:
				f_sel.append(f.index)

		# Store
		mat = context.object.material_slots[ob_seg.active_material_index].material
		mat_seg_dict[mat.name] = (f_sel,mat)

		# Next!
		ob_seg.active_material_index += 1

	# Deselect this object
	ob_seg.select = False

	print("Done")

	############################
	############################
	# Make the compartment objects
	############################
	############################

	# Turn on edit mode
	bpy.ops.object.mode_set(mode='EDIT')

	# Go through all the surface regions
	for reg_surf in reg_list_surf:
		name_surf = reg_surf.name
		if len(name_surf) < 3 or name_surf[:3] != "sc_":
			continue
		if len_surf == 15:
			sc_id = (int(name_surf[3:5]),int(name_surf[6:8]))
		elif len_surf == 21:
			sc_id = (int(name_surf[3:5]),int(name_surf[6:8]),int(name_surf[12:14]))

		print("Making object for compartment: " + str(sc_id))

		# Reset the vertex, edge, face lists
		v_new_list = []
		v_id_surf_dict = {} # From the surface idxs to the new vert list
		v_id_seg_dict = {} # From the segment idxs to the new vert list
		e_new_list = []
		f_new_list = []

		# Dictionary from face indexes in the old object to the new one
		f_id_surf_dict = {}
		f_id_seg_dict = {}

		# Get the faces in this region
		f_reg_surf_list = list(reg_surf.get_region_faces(ob_surf.data))

		# Add to the lists

		# Vertices
		f_v_list = [ob_surf.data.polygons[item].vertices for item in f_reg_surf_list]
		v_unique_list = [item for sublist in f_v_list for item in sublist]
		v_unique_list = list(set(v_unique_list))
		for v in v_unique_list:
			v_new_list.append(tuple(ob_surf.data.vertices[v].co))
			v_id_surf_dict[v] = len(v_new_list) - 1

		# Edges, faces
		for i_fs,fs in enumerate(f_v_list):
			# Face
			nfs = [v_id_surf_dict[v] for v in fs]
			f_new_list.append(nfs)

			# Face index mapping
			f_id_surf_dict[f_reg_surf_list[i_fs]] = len(f_new_list)-1

			# Edges
			e1 = sorted([nfs[0],nfs[1]])
			e2 = sorted([nfs[0],nfs[2]])
			e3 = sorted([nfs[1],nfs[2]])
			for ei in [e1,e2,e3]:
				if not ei in e_new_list:
					e_new_list.append(ei)

		# Go through all segments
		for reg_seg in reg_list_seg:
			if len(reg_seg.name) < 3 or reg_seg.name[:3] != "sc_":
				continue
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
				f_reg_seg_list = list(reg_seg.get_region_faces(ob_seg.data))

				# Add to the lists

				# Vertices
				f_v_list = [ob_seg.data.polygons[item].vertices for item in f_reg_seg_list]
				v_unique_list = [item for sublist in f_v_list for item in sublist]
				v_unique_list = list(set(v_unique_list))
				v_unique_co_list = [tuple(ob_seg.data.vertices[v].co) for v in v_unique_list]
				# Check if it already exists
				for i_v, vco in enumerate(v_unique_co_list):
					v = v_unique_list[i_v]
					if not v in v_id_seg_dict:
						# Search if the vertex already exists from the surface
						try:
							idx = v_new_list.index(vco)
							v_id_seg_dict[v] = idx 
						# Append it as a new vertex
						except ValueError:
							v_new_list.append(vco)
							v_id_seg_dict[v] = len(v_new_list) - 1

				# Edges, faces
				for i_fs,fs in enumerate(f_v_list):
					# Face
					nfs = [v_id_seg_dict[v] for v in fs]
					f_new_list.append(nfs)

					# Face index mapping
					f_id_seg_dict[f_reg_seg_list[i_fs]] = len(f_new_list)-1

					# Edges
					e1 = sorted([nfs[0],nfs[1]])
					e2 = sorted([nfs[0],nfs[2]])
					e3 = sorted([nfs[1],nfs[2]])
					for ei in [e1,e2,e3]:
						if not ei in e_new_list:
							e_new_list.append(ei)

		###################################
		###################################
		# Finally: make the object
		###################################
		###################################

		# Verion 1:
		# Make it using bmesh....
		'''
		bm = bmesh.new()
		v_arr = []
		for v in v_new_list:
			v_arr.append(bm.verts.new(v))
		f_arr = []
		for f in f_new_list:
			bm.faces.new((v_arr[f[0]], v_arr[f[1]], v_arr[f[2]]))

		if len_seg == 19:
			obj_new_name =  "sc_%02d_%02d_compartment" % sc_id
		elif len_seg == 31:
			obj_new_name = "sc_%02d_%02d_sg_%02d_compartment" % sc_id

		mesh_new = bpy.data.meshes.new(obj_new_name  + "_mesh")
		bm.to_mesh(mesh_new)
		bm.free()
		'''

		# Verion 2:
		# Make it using from_pydata....
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

		# Fill holes
		# Where do they come from? Why are they invisible?
		# Nobody fucking knows but apparently this makes everything ok
		bpy.ops.object.mode_set(mode='OBJECT')
		for obj0 in bpy.data.objects: # Deselect all
			obj0.select = False
		obj_new.select = True
		context.scene.objects.active = obj_new
		bpy.ops.object.mode_set(mode='EDIT')
		bpy.ops.mesh.select_all(action='SELECT')
		bpy.ops.mesh.fill_holes()
		bpy.ops.mesh.normals_make_consistent(inside=False) # Also fix the normals
		bpy.ops.object.mode_set(mode='OBJECT')

		###################################
		###################################
		# Copy over ALL MCell regions to the object
		###################################
		###################################	

		# Dictionary of regions to make on the new object
		# Name : face idxs
		new_reg_dict = {}

		# Go through all the surface regions
		for reg_surf in reg_list_surf:
			# Get the faces
			f_surf_list = reg_surf.get_region_faces(ob_surf.data)

			# Find if there are any faces in this original region that also exist on the new object
			# Go through all faces
			for f in f_surf_list:
				if f in f_id_surf_dict:
					# Add to the dict of regions to make
					if reg_surf.name in new_reg_dict:
						new_reg_dict[reg_surf.name].append(f_id_surf_dict[f])
					else:
						new_reg_dict[reg_surf.name] = [f_id_surf_dict[f]]

		# Go through all the segment regions
		for reg_seg in reg_list_seg:
			# Get the faces
			f_seg_list = reg_seg.get_region_faces(ob_seg.data)

			# Find if there are any faces in this original region that also exist on the new object
			# Go through all faces
			for f in f_seg_list:
				if f in f_id_seg_dict:
					# Add to the dict of regions to make
					if reg_seg.name in new_reg_dict:
						new_reg_dict[reg_seg.name].append(f_id_seg_dict[f])
					else:
						new_reg_dict[reg_seg.name] = [f_id_seg_dict[f]]

		# Make the regions

		# Ensure again the new obj is active
		context.scene.objects.active = obj_new

		# Add each of the surfaces as a region
		for reg_name, f_list in new_reg_dict.items():

			# Make region
			obj_new.mcell.regions.add_region_by_name(context,reg_name)

			# Get region
			new_reg = obj_new.mcell.regions.region_list[reg_name]

			# Assign faces
			new_reg.set_region_faces(obj_new.data, f_list)

		# Update (but dont validate because who knows)
		obj_new.data.update()

		###################################
		###################################
		# Add the material to the object
		###################################
		###################################

		# Ensure YET again the new obj is active
		context.scene.objects.active = obj_new

		# Go through all material face lists on both objects
		for dct_tple in [(f_id_surf_dict,mat_surf_dict),(f_id_seg_dict,mat_seg_dict)]:
			f_id_dict = dct_tple[0]
			mat_dict = dct_tple[1]
			for mat_name,tple in mat_dict.items():

				# Get out of the tuple
				mat_face_list = tple[0]
				mat = tple[1]

				# Deselect all vertices
				bpy.ops.object.mode_set(mode='EDIT')
				bpy.ops.mesh.select_all(action='DESELECT')

				# Select all vertices
				ANY_FLAG = False
				for f in mat_face_list:
					if f in f_id_dict:
						obj_new.data.polygons[f_id_dict[f]].select = True
						ANY_FLAG = True

				bpy.ops.object.mode_set(mode='OBJECT')

				# Only proceed if at least some faces of this material are on this object
				if ANY_FLAG == True:

					# Add a material slot
					bpy.ops.object.material_slot_add()

					# Assign a material to the last slot
					context.object.material_slots[context.object.material_slots.__len__() - 1].material = mat

					# Assign the material on the selected vertices
					bpy.ops.object.material_slot_assign()

	print("> Finished: f_regions_to_compartments")







