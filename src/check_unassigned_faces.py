import bpy, bmesh
from mathutils import Vector
import math

# Main

def f_check_unassigned_faces(context, RETURN_FLAG):
	
	if not RETURN_FLAG:
		print("> Running: f_check_unassigned_faces")

	# Get the active object
	ob_list = context.selected_objects
	
	if len(ob_list) == 1:
		ob = ob_list[0]
	
		# Get the region list
		reg_list = ob.mcell.regions.region_list

		# Store all the face indexes
		f_list = []

		# Go through the regions
		for reg in reg_list:
		
			# Get the faces
			f_reg_list = reg.get_region_faces(ob.data)

			# Store
			f_list += f_reg_list

		# Get the list of faces that exist
		f_real_list = list(range(0,len(ob.data.polygons)))

		# Which are missing
		f_missing = list(set(f_real_list)-set(f_list))

		# Empty?
		if len(f_missing) == 0:
			print("No unassigned faces: " + str(len(f_real_list)) + " / " + str(len(f_list)))
		else:
			print("Number of unassigned faces: " + str(len(f_missing)))

			# Show the faces by selecting
			if not RETURN_FLAG:

				# Make sure everything is de-selected before we start
				bpy.ops.object.mode_set(mode='EDIT')
				bpy.ops.mesh.select_all(action='DESELECT')
				bpy.ops.object.mode_set(mode='OBJECT')

				for f in f_missing:
					ob.data.polygons[f].select = True

				# Return in edit mode to show
				bpy.ops.object.mode_set(mode='EDIT')

			else:

				# Return flag is true, so return the list of missing faces
				return f_missing

	if not RETURN_FLAG:
		print("> Finished: f_check_unassigned_faces")

def f_select_unassigned_linked_faces(context):

	print("> Running: f_select_unassigned_linked_faces")

	# Get the active object
	ob_list = context.selected_objects
	
	if len(ob_list) == 1:
		ob = ob_list[0]

		# Check that only one face is currently selected
		if ob.data.total_face_sel != 1:
			raise TypeError("Please select only one face.")

		# Get the current face
		bpy.ops.object.mode_set(mode='OBJECT')
		face_index = [f.index for f in ob.data.polygons if f.select][0]

		# Get the list of unassigned faces
		f_missing = f_check_unassigned_faces(context, True)

		# Check that the selected face is in here
		if not face_index in f_missing:
			raise TypeError("Selected face is not unassigned.")

		# Get the edges of all these missing faces
		e_missing = [set(list(ob.data.polygons[i].edge_keys)) for i in f_missing]

		# Go through and find linked faces
		f_check = [face_index]
		f_done = []

		while len(f_check) > 0:
			# Face to check
			face_index = f_check[0]

			# Get edges
			edges = set(list(ob.data.polygons[face_index].edge_keys))

			# Check against other edges
			for i_f, edges_missing in enumerate(e_missing):
				if len(edges - edges_missing) != 3:
					# An edge is shared => linked face
					f_linked = f_missing[i_f]
					# Check that it hasn't already been done and isn't myself
					if not f_linked in f_done and not face_index == f_linked:
						# Select
						ob.data.polygons[f_linked].select = True
						# Store for checking
						f_check.append(f_linked)

			# Make sure not to double count
			f_done.append(face_index)

			# Delete
			del f_check[0]

		# Return in edit mode
		bpy.ops.object.mode_set(mode='EDIT')

	print("> Finished: f_select_unassigned_linked_faces")
