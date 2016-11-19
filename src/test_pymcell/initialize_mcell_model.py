# import bpy, bmesh
# from mathutils import Vector
import math

import sys

# For CellBlender
# import cellblender

# For pymcell
import pymcell as m

# For importing MDL geometries
import import_mdl_geometry

# For computing necessary statistics about the geometry
import calculate_geometry_info

# Multiprocessing to run MCell commands (else: Blender crashes)
from multiprocessing import Process

# For comparing names
import re

# Class to hold MCell molecule properties AND the actual MCell object
class MCell_Molecule:

	def __init__(self, *args, **kwargs):

		# Molecule
		if 'name' in kwargs:
			self.name = kwargs['name']
		else:
			self.name = ""

		# Type
		if 'type' in kwargs:
			self.type = kwargs['type']
		else:
			self.type = ""

		# Object name
		if 'geom_name' in kwargs:
			self.geom_name = kwargs['geom_name']
		else:
			geom_name = ""

		# Region id
		if 'reg_id' in kwargs:
			self.reg_id = kwargs['reg_id']
		else:
			self.reg_id = ()

		# Diffusion constant
		if 'dc' in kwargs:
			self.dc = kwargs['dc']
		else:
			self.dc = 0

		# Number to release
		if 'nrel' in kwargs:
			self.nrel = kwargs['nrel']
		else:
			self.nrel = 0

		# Surface flag
		if 'surf_flag' in kwargs:
			self.surf_flag = kwargs['surf_flag']
		else:
			self.surf_flag = False

		self.mcell_obj = None


# Class to hold geometry info
class MCell_Geometry:

	def __init__(self, *args, **kwargs):

		# Object name
		if 'name' in kwargs:
			self.name = kwargs['name']
		else:
			self.name = ""

		# List of vertices
		if 'seg_vert_list' in kwargs:
			self.seg_vert_list = kwargs['seg_vert_list']
		else:
			self.seg_vert_list = []
		if 'surf_vert_list' in kwargs:
			self.surf_vert_list = kwargs['surf_vert_list']
		else:
			self.surf_vert_list = []
		if 'exc_vert_list' in kwargs:
			self.exc_vert_list = kwargs['exc_vert_list']
		else:
			self.exc_vert_list = []

		# List of elements
		if 'seg_el_list' in kwargs:
			self.seg_el_list = kwargs['seg_el_list']
		else:
			self.seg_el_list = []
		if 'surf_el_list' in kwargs:
			self.surf_el_list = kwargs['surf_el_list']
		else:
			self.surf_el_list = []
		if 'exc_el_list' in kwargs:
			self.exc_el_list = kwargs['exc_el_list']
		else:
			self.exc_el_list = []

		# Dict of regions
		if 'seg_reg_dict' in kwargs:
			self.seg_reg_dict = kwargs['seg_reg_dict']
		else:
			self.seg_reg_dict = {}
		if 'surf_reg_dict' in kwargs:
			self.surf_reg_dict = kwargs['surf_reg_dict']
		else:
			self.surf_reg_dict = {}
		if 'exc_reg_dict' in kwargs:
			self.exc_reg_dict = kwargs['exc_reg_dict']
		else:
			self.exc_reg_dict = {}

		# Dict of number of segments per section
		if 'nsg_sc_dict' in kwargs:
			self.nsg_sc_dict = kwargs['nsg_sc_dict']
		else:
			self.nsg_sc_dict = {}

		# Dict of surface ares
		if 'surf_area_dict' in kwargs:
			self.surf_area_dict = kwargs['surf_area_dict']
		else:
			self.surf_area_dict = {}

		# The mcell objects
		self.seg_obj = None
		self.surf_obj = None
		self.exc_obj = None


# Class to hold a release site
class MCell_Rel_Site:

	def __init__(self, *args, **kwargs):

		# Release site (object) name
		if 'name' in kwargs:
			self.name = kwargs['name']
		else:
			self.name = ""

		# Shape
		if 'shape' in kwargs:
			self.shape = kwargs['shape']
		else:
			self.shape = ""

		# Object name
		if 'obj_name' in kwargs:
			self.obj_name = kwargs['obj_name']
		else:
			self.obj_name = ""

		# Location
		if 'location' in kwargs:
			self.location = kwargs['location']
		else:
			self.location = ()

		# Molecule Type
		if 'mol_type' in kwargs:
			self.mol_type = kwargs['mol_type']
		else:
			self.mol_type = ""

		# Number to release
		if 'nrel' in kwargs:
			self.nrel = kwargs['nrel']
		else:
			self.nrel = 0

		# Site diameter
		if 'diam' in kwargs:
			self.diam = kwargs['diam']
		else:
			self.diam = 0

		# Release delay
		if 'delay' in kwargs:
			self.delay = kwargs['delay']
		else:
			self.delay = 0

		# Release object
		self.rel_obj = None

###
# Function to initialize the MCell
###

def tmp():
	import torus

	world = m.mcell_create()
	m.mcell_init_state(world)

	dt = 1e-5
	iterations = 100
	m.mcell_set_time_step(world, dt)
	m.mcell_set_iterations(world, iterations)

	# Define one surface molecule and three volume molecules
	sm1_sym = m.create_species(world, "sm1", 1e-6, True)
	vm1_sym = m.create_species(world, "vm1", 1e-6, False)
	vm2_sym = m.create_species(world, "vm2", 1e-6, False)
	vm3_sym = m.create_species(world, "vm3", 1e-6, False)

	# Define reactions
	# vm1 + vm2 -> vm3 [1e8]
	reactants1 = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
	reactants1 = m.mcell_add_to_species_list(vm2_sym, False, 0, reactants1)
	products1 = m.mcell_add_to_species_list(vm3_sym, False, 0, None)
	m.create_reaction(world, reactants1, products1, 1e8)
	# vm3 -> NULL [1e5]
	reactants2 = m.mcell_add_to_species_list(vm3_sym, False, 0, None)
	m.create_reaction(world, reactants2, None, 0.01, name="rxn")

	scene_name = "Scene"
	scene = m.create_instance_object(world, scene_name)

	# Create a spherical release site
	pos_vec3 = m.Vector3()
	diam_vec3 = m.Vector3(0.015, 0.015, 0.015)
	# XXX: It seems to be necessary to return some or all of these objects in
	# order to have a functioning release site even though we don't use them
	# anywhere after this call.
	position, diameter, sphere_release_object = m.create_release_site(
		world, scene, pos_vec3, diam_vec3, m.SHAPE_SPHERICAL, 500, vm1_sym,
		"vm1_rel")
	pos_vec3b = m.Vector3(0.05, 0.05, 0.00)
	diam_vec3b = m.Vector3(0.025, 0.025, 0.05)
	position2, diameter2, cube_release_object = m.create_release_site(
		world, scene, pos_vec3b, diam_vec3b, m.SHAPE_CUBIC, 500, vm2_sym,
		"vm2_rel")

	# Create box object
	box_name = "Box"
	box_mesh = m.create_box(world, scene, 0.1, box_name)

	box_region_name = "side"
	box_face_list = [1,2]
	box_region = m.create_surface_region(
		world, box_mesh, box_face_list, box_region_name)

	torus_name = "Torus"
	torus_mesh = m.create_polygon_object(
		world, torus.vert_list, torus.face_list, scene, torus_name)

	# Create surface region on half the torus
	# XXX: Creating a region is currently required when creating torus_mesh objects
	torus_region_name = "half_torus"
	torus_region = m.create_surface_region(
		world, torus_mesh, torus.surf_reg_face_list, torus_region_name)

	region_release_object = m.create_region_release_site(
		world, scene, torus_mesh, "vm1_torus_rel", "ALL", 1000, vm1_sym)

	# create surface class
	sc_sm1_sym = m.create_surf_class(world, "sc_release_y")
	# create releases using a surface class (i.e. not a release object)
	# mdl equivalent: MOLECULE_DENSITY {sm1' = 1000}
	sm1 = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
	smd = m.mcell_add_mol_release_to_surf_class(
		world, sc_sm1_sym, sm1, 1000, 0, None)

	# create surface class that is reflective to sm1
	m.mcell_add_surf_class_properties(world, m.RFLCT, sc_sm1_sym, sm1_sym, 0)
	# m.mcell_add_surf_class_properties(world, m.SINK, sc_sm1, sm1_sym, 0)
	m.mcell_assign_surf_class_to_region(sc_sm1_sym, torus_region)

	m.mcell_delete_species_list(sm1)

	# Create reaction on sc_sm1 sm1, -> sm1'
	sc_surf = m.mcell_add_to_species_list(sc_sm1_sym, True, 1, None)
	reactantsSurf = m.mcell_add_to_species_list(sm1_sym, True, 1, None)
	productsSurf = m.mcell_add_to_species_list(sm1_sym, True, -1, None)
	m.create_reaction(
		world, reactantsSurf, productsSurf, 1e4, surf_class=sc_surf, name="rxnSurf")

	# Create viz data
	viz_list = m.mcell_add_to_species_list(vm1_sym, False, 0, None)
	viz_list = m.mcell_add_to_species_list(vm2_sym, False, 0, viz_list)
	viz_list = m.mcell_add_to_species_list(vm3_sym, False, 0, viz_list)
	viz_list = m.mcell_add_to_species_list(sm1_sym, True, 0, viz_list)
	m.mcell_create_viz_output(
		world, "./viz_data/test", viz_list, 0, iterations, 1)

	# Create reaction data
	box_sym = m.mcell_get_obj_sym(box_mesh)
	count_list1, os1, out_times1, output1 = m.create_count(
		world, box_sym, vm1_sym, "react_data/vm1_%s.dat" % box_name, 1e-5)
	torus_reg_sym = m.mcell_get_reg_sym(torus_region)
	count_list2, os2, out_times2, output2 = m.create_count(
		world, torus_reg_sym, sm1_sym, "react_data/sm1_reg.dat", dt)
	count_list3, os3, out_times3, output3 = m.create_count(
		world, None, vm1_sym, "react_data/vm1_world.dat", dt)
	count_list4, os4, out_times4, output4 = m.create_count(
		world, box_sym, vm3_sym, "react_data/vm3_%s.dat" % box_name, 1e-5)

	m.mcell_init_simulation(world)
	m.mcell_init_output(world)

	output_freq = 10
	for i in range(iterations):
		vm3_count = m.mcell_get_count(
			"vm3", "%s.%s,ALL" % (scene_name, box_name), world)
		# print(vm3_count)
		# When vm3 hits some arbitrary threshold value (i.e. 400), ramp up the
		# rate constant of vm3->NULL. This is just a simple test, but we'll
		# need to do something analagous when interfacing with pyNEURON.
		if (vm3_count > 400):
			m.mcell_modify_rate_constant(world, "rxn", 1e8)
		m.mcell_run_iteration(world, output_freq, 0)

	m.mcell_flush_data(world)
	m.mcell_print_final_warnings(world)
	m.mcell_print_final_statistics(world)


def simple():
	print("hello")

def f_init_mcell(DIR, timestep_info, MCELL_geom_list, MCELL_mol_obj_dict, MCELL_rel_site_list):

	# Unbundle the timestep info
	dt_M, n_iter_N, n_steps_M = timestep_info

	# Purge MCELL_mol_obj_dict to have only molecules by type
	mol_types_done = []
	MCELL_mol_obj_types_dict = {}
	for mol in MCELL_mol_obj_dict.values():
		if mol.type in mol_types_done:
			continue

		# Add
		MCELL_mol_obj_types_dict[mol.type] = mol

		# Fix name
		MCELL_mol_obj_types_dict[mol.type].name = mol.type

		# Don't repeat this type of molecule
		mol_types_done.append(mol.type)

	# Create the world
	world_main = m.mcell_create()
	m.mcell_init_state(world_main)

	# Timestep, iters
	m.mcell_set_time_step(world_main, dt_M)
	m.mcell_set_iterations(world_main, 1)

	###
	# Define molecules by TYPE, i.e. not for each segment
	###
	'''
	for mol in MCELL_mol_obj_types_dict.values():

		# Make the molecule
		tmp_particle = m.create_species(world_main, mol.type, mol.dc, mol.surf_flag)
		mol.mcell_obj = tmp_particle

	###
	# Create a scene
	###

	scene_name = "Scene"
	scene = m.create_instance_object(world_main, scene_name)

	###
	# Create geometry
	###

	for geom in MCELL_geom_list:
		# Dendrites
		if geom.name[0] == 'd':
			geom.exc_obj = m.create_polygon_object(world_main, geom.exc_vert_list, geom.exc_el_list, scene, geom.name + "_Excluded")

	###
	# Create release sites
	###

	# The mesh to release into
	rel_geom_obj = [geom.exc_obj for geom in MCELL_geom_list if geom.name[0] == 'd'][0]
	rel_geom_sym = m.mcell_get_obj_sym(rel_geom_obj)

	# Something to store the counts in
	store_count_dict = {}

	for rel in MCELL_rel_site_list:
		# Only release molecules that are meant to be released initially
		if rel.delay == 0:
			# The molecule object
			rel_mol = [mol for mol in MCELL_mol_obj_types_dict.values() if mol.type == rel.mol_type][0]
			rel.rel_obj = m.create_region_release_site(world_main, scene, rel_geom_obj, rel.name, "ALL", rel.nrel, rel_mol.mcell_obj)

			###
			# Create output data
			###

			store_count_dict[rel_mol.name] = m.create_count(world_main, rel_geom_sym, rel_mol.mcell_obj, DIR + "react_data/"+rel_mol.name+".dat", dt_M)

	'''

	# Init
	m.mcell_init_simulation(world_main)
	m.mcell_init_output(world_main)

	# Run
	output_freq = 1 # Output at every timestep, rather than at every "output_freq"-th
	m.mcell_run_iteration(world_main, output_freq, 1)

	# Finish
	m.mcell_flush_data(world_main)
	m.mcell_print_final_warnings(world_main)
	m.mcell_print_final_statistics(world_main)



###
# Main
###


if __name__ == "__main__":

	print('> Running: initialize_mcell_model.py')

	'''

	# Working dir
	global DIR
	DIR = "/home/oernst/Research/cnl/neuron_mcell/neuron_launcher/src/test_pymcell/"

	##########
	# Read geometry data
	##########


	# List to store geometries
	MCELL_geom_list = []

	# List of objects to import
	obj_name_list = ['a170','d001']

	# Import objects
	for obj_name in obj_name_list:

		# Surface
		surf_name = obj_name + "_Surface"
		mn_surface_vert_list, mn_surface_el_list, mn_surface_reg_dict = import_mdl_geometry.import_mdl(DIR+surf_name+".mdl", surf_name)

		# Segment
		seg_name = obj_name + "_Segment"
		mn_segment_vert_list, mn_segment_el_list, mn_segment_reg_dict = import_mdl_geometry.import_mdl(DIR+seg_name+".mdl", seg_name)

		# If it is a dendrite, we need to exclude some regions of the volume
		if obj_name[0] == 'd':

			# Exclude
			exc_name = obj_name + "_Excluded"
			mn_exclude_vert_list, mn_exclude_el_list, mn_exclude_reg_dict = import_mdl_geometry.import_mdl(DIR+exc_name+".mdl", exc_name)

		# Get the number of segments per section
		mn_nsg_sc_dict = calculate_geometry_info.compute_nsg_sc(mn_surface_reg_dict)

		# Get the surface areas
		mn_surf_area_dict = calculate_geometry_info.compute_surface_areas(mn_surface_vert_list, mn_surface_el_list, mn_surface_reg_dict)

		# Make geometry object
		tmp_geom = MCell_Geometry(
			name=obj_name, 
			surf_vert_list=mn_surface_vert_list,
			surf_el_list=mn_surface_el_list,
			surf_reg_dict=mn_surface_reg_dict,
			seg_vert_list=mn_segment_vert_list,
			seg_el_list=mn_segment_el_list,
			seg_reg_dict=mn_segment_reg_dict,
			nsg_sc_dict=mn_nsg_sc_dict,
			surf_area_dict=mn_surf_area_dict
			)

		# If it is a dendrite, set the excluded volume
		if obj_name[0] == 'd':
			tmp_geom.exc_vert_list=mn_exclude_vert_list,
			tmp_geom.exc_el_list=mn_exclude_el_list,
			tmp_geom.exc_reg_dict=mn_exclude_reg_dict,

		# Add to list
		MCELL_geom_list.append(tmp_geom)


	##########
	# NEURON: Timestep/iterations
	##########


	# Timestep
	dt_N = 0.025 # [ms]
	
	# Number iterations
	n_iter_N = 200


	##########
	# MCell: Timestep/iterations
	##########	


	# Number of MCell steps to take for each NEURON timestep
	n_steps_M = 10
	dt_M = dt_N / n_steps_M
	# m.mcell_set_time_step(world_main, dt_M)
	
	# Number iterations
	# m.mcell_set_iterations(world_main, n_iter_N * n_steps_M)

	##########
	# MCell: Read parameters from CellBlender into a dict
	##########

	# Read into dict
	cb_param_dict = {}
	for this_param_dict in bpy.context.scene.mcell.parameter_system['gp_dict'].values():
		cb_param_dict[this_param_dict['name']] = this_param_dict['value']


	##########
	# MCell: Molecules: discover information
	##########


	# Store all the molecules
	MCELL_mol_obj_dict = {}

	# Get the list of molecules in CellBlender
	cb_mol_list = bpy.context.scene.mcell.molecules.molecule_list

	# Make all the molecule objects for each region and store them
	# Go through all molecules
	for mol in cb_mol_list:

		###
		# Exclude the following molecules from copying
		###

		# Skip glutamate - it is extracellular
		if mol.name == "Glu":
			continue

		###
		# Make copies of the following mols
		###

		# print("Making copies of molecule: " + str(mol.name))

		# Read the mol info

		mol_dc = mol.diffusion_constant.get_expr()
		try:
			mol_dc = float(mol_dc)
		except:
			if mol_dc in cb_param_dict:
				# print("Assigning val: " + str(cb_param_dict[mol_dc]) + " to dc: " + mol_dc)
				mol_dc = cb_param_dict[mol_dc]
			else:
				print("Error: could not find definition for DC: " + str(mol_dc))
				mol_dc = 0.0

		mol_vol_surf_type = str(mol.type)
		if mol_vol_surf_type == '3D':
			mol_surf_flag = False
		else:
			mol_surf_flag = True

		# Make a copy of this molecule for all segments for all sections for all objects....
		# Go through all objects
		for geom in MCELL_geom_list:
			for sc_id, nsg in geom.nsg_sc_dict.items():
				for isg in range(1,nsg+1):

					mol_name = mol.name + "_" + geom.name + "_sc_%02d_%02d_sg_%02d"%(sc_id[0],sc_id[1],isg)

					# Make the mol
					tmp_obj = MCell_Molecule(
						name = mol_name, 
						type = mol, 
						geom_name = geom.name,
						reg_id = (sc_id[0],sc_id[1],isg), 
						dc = mol_dc, 
						surf_flag = mol_surf_flag)

					# Store
					MCELL_mol_obj_dict[mol_name] = tmp_obj


	##########
	# MCell: Discover release sites
	##########


	# Store all release sites
	MCELL_rel_site_list = []


	###
	# Release onto the particular synapse for Glutamate
	###

	# Release site name pattern
	rel_name_pattern = 'd001p_sy.*_rs'

	# Go through all objects
	for ob in bpy.data.objects:

		# Check the name
		if re.match(rel_name_pattern, ob.name) != None:

			# This is a valid release site
			print("Found a valid release site: " + str(ob.name))

			# Location
			rel_loc = tuple(ob.data.vertices[0].co)

			# Number to release
			rel_nrel = cb_param_dict["n_Glu"]

			# Release delay
			rel_delay = cb_param_dict["glu_release_delay"]

			# Make a new release site
			tmp_rel = MCell_Rel_Site(
				name=ob.name,
				shape="SPHERICAL",
				location=rel_loc,
				mol_type="Glu",
				nrel=rel_nrel,
				diam=0,
				delay=rel_delay
				)

			# Store the release site information
			MCELL_rel_site_list.append(tmp_rel)

	###
	# Release into the excluded volume
	###

	# Tuples of info
	rel_info = [('SPINECa', 'Ca', 'basalCa'),
		('SPINEFASTB','fast_sp.B','fast_cbp_concentration * bound_cbp_fraction_spine'),
		('SPINEFASTU','fast_sp.U', 'fast_cbp_concentration * unbound_cbp_fraction_spine'),
		('SPINEMEDB', 'medium_sp.B', 'medium_cbp_concentration * bound_cbp_fraction_spine'),
		('SPINEMEDU', 'medium_sp.U','medium_cbp_concentration * unbound_cbp_fraction_spine'),
		('SPINESLOWB','slow_sp.B', 'slow_cbp_concentration * bound_cbp_fraction_spine'),
		('SPINESLOWU', 'slow_sp.U', 'slow_cbp_concentration * unbound_cbp_fraction_spine'),
		('fluo4_unbound_rel', 'fluo4.U', 'fluo4_conc * fluo4_unbound_fraction'),
		('fluo4_bound_rel', 'fluo4.B', 'fluo4_conc * fluo4_bound_fraction'),
		('ogb1_unbound_rel', 'ogb1.U', 'ogb1_conc * ogb1_unbound_fraction'),
		('ogb1_bound_rel', 'ogb1.B', 'ogb1_conc * ogb1_bound_fraction'),
		('SPINECALBH0M0', 'calbindin.high0medium0', 'calbindin_concentration * calbindin_H0M0_fraction'),
		('SPINECALBH0M1', 'calbindin.high0medium1', 'calbindin_concentration * calbindin_H0M1_fraction'),
		('SPINECALBH0M2', 'calbindin.high0medium2', 'calbindin_concentration * calbindin_H0M2_fraction'),
		('SPINECALBH1M0', 'calbindin.high1medium0', 'calbindin_concentration * calbindin_H1M0_fraction'),
		('SPINECALBH1M1', 'calbindin.high1medium1', 'calbindin_concentration * calbindin_H1M1_fraction'),
		('SPINECALBH1M2', 'calbindin.high1medium2', 'calbindin_concentration * calbindin_H1M2_fraction'),
		('SPINECALBH2M0', 'calbindin.high2medium0', 'calbindin_concentration * calbindin_H2M0_fraction'),
		('SPINECALBH2M1', 'calbindin.high2medium1', 'calbindin_concentration * calbindin_H2M1_fraction'),
		('SPINECALBH2M2', 'calbindin.high2medium2', 'calbindin_concentration * calbindin_H2M2_fraction')]

	# Release object name
	rel_obj_name = "d001_Exclude"

	# Make the release objects
	for rel in rel_info:

		# Release site name
		rel_site_name = rel[0]

		# Molecule
		rel_mol = rel[1]

		# Value
		rel_val_str = rel[2]
		rel_val_split = rel_val_str.split()
		if len(rel_val_split) == 1:
			rel_val = cb_param_dict[rel_val_str]
		elif len(rel_val_split) == 3:
			rel_val = cb_param_dict[rel_val_split[0]] * cb_param_dict[rel_val_split[2]]
		else:
			print("Error: could not determine the number to release for release site: " + str(rel_site_name))
			continue

		# Make a new release site
		tmp_rel = MCell_Rel_Site(
			name=rel_site_name,
			shape="REGION",
			obj_name=rel_obj_name,
			mol_type=rel_mol,
			nrel=rel_val
			)

		# Store
		MCELL_rel_site_list.append(tmp_rel)

	'''

	##########
	# MCell: Run MCell
	##########

	print("> Running MCell")

	# Bundle the timestep info
	# timestep_info = (dt_M, n_iter_N, n_steps_M)

	p = Process(target=tmp, args=())
	# p = Process(target=f_init_mcell, args=(DIR, timestep_info, MCELL_geom_list, MCELL_mol_obj_dict, MCELL_rel_site_list,))
	p.start()
	p.join()

	print('> Finished: initialize_mcell_model.py')


