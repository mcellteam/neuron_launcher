import os

from math import *

import sys

import random

import copy

# For pymcell
import pymcell as m

# For NEURON
from neuron import h

# For importing MDL geometries
import import_mdl_geometry

# For computing necessary statistics about the geometry
import calculate_geometry_info

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

# Class to hold MCell reactions
class MCell_Reaction:

	def __init__(self, *args, **kwargs):

		# Name
		if 'name' in kwargs:
			self.name = kwargs['name']
		else:
			self.name = ""
		
		# RC forward
		if 'rc_f' in kwargs:
			self.rc_f = kwargs['rc_f']
		else:
			self.rc_f = 0

		# RC backward
		if 'rc_b' in kwargs:
			self.rc_b = kwargs['rc_b']
		else:
			self.rc_b = 0

		# RC functions
		if 'f_rc_f' in kwargs:
			self.f_rc_f = kwargs['f_rc_f']
		else:
			self.f_rc_f = None
		if 'f_rc_b' in kwargs:
			self.f_rc_b = kwargs['f_rc_b']
		else:
			self.f_rc_b = None

		# Reaction class, i.e. real or compartment
		'''
		if 'rxn_type' in kwargs:
			self.rxn_type = kwargs['rxn_type']
		else:
			self.rxn_type = ""
		'''

		# Types of molecules, i.e. name of the molecules (WITHOUT location id) on either side
		if 'mol_types' in kwargs:
			self.mol_types = kwargs['mol_types']
		else:
			self.mol_types = []

		# Region name in which this reaction occurs
		if 'reg_id' in kwargs:
			self.reg_id = kwargs['reg_id']
		else:
			self.reg_id = ""

		# Reactant list
		self.reactant_list = None

		# Product list
		self.product_list = None

		# Surface list
		self.surf_list = None

# Read an SWC file, make a NEURON model
def neuron_model_from_swc(swc_filepath, nsg_sc_dict):
	
	# Store the neuron sections
	nrn_sc_dict = {}

	# Store the neuron segments created
	nrn_sg_dict = {}

	# Read the SWC file
	swc_data = []
	f = open(swc_filepath,'r')
	for i_line,line in enumerate(f):
		line_s = line.split();
		if len(line_s) > 0 and line_s[0] != "#":
			line_data = []
			line_data.append(int(float(line_s[0])))
			line_data.append(int(float(line_s[1])))
			line_data.append(float(line_s[2]))
			line_data.append(float(line_s[3]))
			line_data.append(float(line_s[4]))
			line_data.append(float(line_s[5])) # Radius
			line_data.append(int(float(line_s[6])))
			swc_data.append(line_data)

	f.close()

	# Find connections
	pt_connect = []
	for i in range(0,len(swc_data)):
		pt_connect.append([])

	for i,data in enumerate(swc_data):
		if i > 0: # The first is the -1 index
			pt_connect[data[6]-1].append(i+1)

	# Make the model
	for i_pt,conn_pts in enumerate(pt_connect):
		pt = i_pt + 1
		
		for i_conn_pt,conn_pt in enumerate(conn_pts):
			min_pt = min(pt,conn_pt)
			max_pt = max(pt,conn_pt)
			min_dat = swc_data[min_pt-1]
			max_dat = swc_data[max_pt-1]
			
			# Section
			sc = h.Section(name="sc_%02d_%02d"%(min_pt,max_pt))
			# print("Creating section: " + str("sc_%02d_%02d"%(min_pt,max_pt)))

			# Section data
			# sc.nseg = 10 # read this out of the mdl file
			sc.L = sqrt(sum([(min_dat[j]-max_dat[j])**2 for j in [2,3,4]])) # um
			sc.cm = 1 # [uF/cm2]

			# Inter-compartment resistivity
			sc.Ra = 20000 # [ohm-cm]

			# Number of segments
			sc.nseg = nsg_sc_dict[(min_pt,max_pt)]
			
			# Diameter
			min_diam = 2.0*min_dat[5]
			max_diam = 2.0*max_dat[5]
			delta_diam = max_diam - min_diam
			for sg in sc:
				sg.diam = min_diam + sg.x*delta_diam # [um]
		
			# Insert mechanisms
			sc.insert("MCell_na")
			sc.insert("MCell_other")

			# Set location inside each mod mechanism - the mod files need to know where they are
			loc_seg = 1;
			for sg in sc:
				# Works because segments are iterated over from x=0 to x=1 ie min to max
				for mech in sg: # iterates over the segment mechanisms
					
					if mech.name() == "MCell_na":
						# Set location
						mech.loc_sec_min = min_pt
						mech.loc_sec_max = max_pt
						mech.loc_seg = loc_seg
						# Set the initial fraction of open channels to the equilibrium values
						mech.frac_open = 8.839737624582281e-05

				# Store the segment
				nrn_sg_dict[(min_pt,max_pt,loc_seg)] = sg

				loc_seg += 1

			# Section connection
			if not (i_pt == 0 and i_conn_pt == 0): # If there are no sections yet, no connections
				# Search for a connecting point
				for other_sc_id in nrn_sc_dict.keys():
					if other_sc_id[0] == pt:
						sc.connect(nrn_sc_dict[other_sc_id](0))
						break
					elif other_sc_id[1] == pt:
						sc.connect(nrn_sc_dict[other_sc_id](1))
						break
		
			# Store the section
			nrn_sc_dict[(min_pt,max_pt)] = sc

	return nrn_sc_dict, nrn_sg_dict

# Generate the reaction list
def generate_rxn_list(mcell_rxn_real_dict, mcell_rxn_comp_list, mol_obj_dict, rxn_list, mol_names, v0, nsg_sc_dict, mn_segment_reg_dict, mn_segment_surf_class_dict):

	print("> Making actual reactions")

	break_flag = False
	# Fix all the reactions
	for rxn in rxn_list:

		# All the sections, segments
		for sec_id,nseg in nsg_sc_dict.items():
			for i_seg in range(1,nseg+1):

				# Region
				reg_id = (sec_id[0],sec_id[1],i_seg)
				reg_name = "sc_%02d_%02d_sg_%02d"%reg_id

				# Mol types
				mol_type_1 = rxn[0]
				mol_type_2 = rxn[3]

				# If backward reaction, flip the types
				if rxn[2] == "<-":
					mol_type_tmp = mol_type_1
					mol_type_1 = mol_type_2
					mol_type_2 = mol_type_tmp

				# Names
				mol_name_1 = mol_type_1 + "_" + reg_name
				mol_name_2 = mol_type_2 + "_" + reg_name
				rxn_name = "rxn_" + mol_name_1 + "_R_" + mol_name_2

				# print("Making reaction: " + str(rxn_name))

				# Reactant/product lists
				reactant_list = m.mcell_add_to_species_list(mol_obj_dict[mol_name_1].mcell_obj, mol_obj_dict[mol_name_1].surf_flag, 0, None)
				product_list = m.mcell_add_to_species_list(mol_obj_dict[mol_name_2].mcell_obj, mol_obj_dict[mol_name_2].surf_flag, 0, None)

				# Make a new rxn object
				# Note that we evaluate the rxn rate at the initial voltage value
				if rxn[2] == "->" or rxn[2] == "<-":
					rxn_obj = MCell_Reaction(name = rxn_name, rc_f = rxn[-1](v0), mol_types = [mol_type_1, mol_type_2], reg_id = reg_id, f_rc_f = rxn[-1])
				elif rxn[2] == "<-->":
					rxn_obj = MCell_Reaction(name = rxn_name, rc_f = rxn[-2](v0), rc_b = rxn[-1](v0), mol_types = [mol_type_1, mol_type_2], reg_id = reg_id, f_rc_f = rxn[-2], f_rc_b = rxn[-1])

				# Set the reactant, product lists
				rxn_obj.reactant_list = reactant_list
				rxn_obj.product_list = product_list

				# Append
				if not reg_id in mcell_rxn_real_dict:
					mcell_rxn_real_dict[reg_id] = [rxn_obj]
				else:
					mcell_rxn_real_dict[reg_id].append(rxn_obj)

	print("> Making inter-compartment reactions")

	# Now we need to do compartment reactions....
	# I.E. Change molecule types when they hit surface boundaries
	# Go through each segment boundary - each segment needs two reactions for EVERY molecule!
	for brdr_name in mn_segment_reg_dict.keys():
		sg_name_bottom = brdr_name[0:14]
		sg_name_top = brdr_name[17:]

		# Go through all of the molecule names
		for mol_name in mol_names:

			mol_name_top = mol_name + "_" + sg_name_top
			mol_name_bottom = mol_name + "_" + sg_name_bottom
			surf_class_name = "SurfClass_"+brdr_name
			rxn_name = "rxn_" + mol_name_top + "_R_" + surf_class_name

			# print("Making reaction: " + str(rxn_name))

			# Reactant/product lists
			reactant_list = m.mcell_add_to_species_list(mol_obj_dict[mol_name_top].mcell_obj, mol_obj_dict[mol_name_top].surf_flag, 0, None)
			product_list = m.mcell_add_to_species_list(mol_obj_dict[mol_name_bottom].mcell_obj, mol_obj_dict[mol_name_bottom].surf_flag, 0, None)
			surf_list = m.mcell_add_to_species_list(mn_segment_surf_class_dict[surf_class_name], True, 1, None)

			# Make 2 rxns IN ONE STATEMENT! :)
			rxn_obj = MCell_Reaction(name = rxn_name, rc_f = 1e50, rc_b = 1e50, mol_types = [mol_name, mol_name])

			# Set the reactant, product lists
			rxn_obj.reactant_list = reactant_list
			rxn_obj.product_list = product_list
			rxn_obj.surf_list = surf_list

			# Append
			mcell_rxn_comp_list.append(rxn_obj)

	return

# Voltage dependent rate constants
def f_alpham(v):
	return 1000.0 * 0.1 * (v + 45.0) / (1.0 - exp(-(v+45.0)/10.0))
def f_betam(v):
	return 1000.0 * 4.0 * exp(-(v+70.0)/18.0);
def f_alphah(v):
	return 1000.0 * 0.07 * exp(-(v+70.0)/20.0);
def f_betah(v):
	return 1000.0 * 1.0 / (1.0 + exp(-(v+40.0)/10.0));
def f_alphan(v):
	return 1000.0 * 0.01 * (v + 60.0) / (1.0 - exp(-(v+60.0)/10.0));
def f_betan(v):
	return 1000.0 * 0.125 * exp(-(v+70.0)/80.0);

# Main

if __name__ == "__main__":

	# Working dir
	global DIR
	DIR = "/home/oernst/Research/cnl/neuron_mcell/neuropil/d001_hh/"

	# Temperature
	celsius = 6.3
		
	##########
	# Read some data
	##########

	# Read the geometry file
	mn_surface_vert_list, mn_surface_el_list, mn_surface_reg_dict = import_mdl_geometry.import_mdl(DIR+ "MN_Surface.mdl", "MN_Surface")
	mn_segment_vert_list, mn_segment_el_list, mn_segment_reg_dict = import_mdl_geometry.import_mdl(DIR+ "MN_Segment.mdl", "MN_Segment")
	
	# Get the number of segments per section
	nsg_sc_dict = calculate_geometry_info.compute_nsg_sc(mn_surface_reg_dict)

	# Get the surface areas
	surf_area_dict = calculate_geometry_info.compute_surface_areas(mn_surface_vert_list, mn_surface_el_list, mn_surface_reg_dict)

	##########
	# NEURON: Make the geometry model
	##########

	nrn_sc_dict, nrn_sg_dict = neuron_model_from_swc(DIR + "d001_v2_cabled.swc", nsg_sc_dict)

	h.topology()
	
	##########
	# MCell: Initialize the MCell world
	##########
	
	world_main = m.mcell_create()
	m.mcell_init_state(world_main)

	##########
	# NEURON: Timestep/iterations
	##########

	# Timestep
	dt_N = 0.025 # [ms]
	h.dt = dt_N
	
	# Number iterations
	n_iter_N = 200

	##########
	# MCell: Timestep/iterations
	##########	

	# Number of MCell steps to take for each NEURON timestep
	n_steps_M = 10
	dt_M = dt_N / n_steps_M
	m.mcell_set_time_step(world_main, dt_M)
	
	# Number iterations
	m.mcell_set_iterations(world_main, n_iter_N * n_steps_M)

	##########
	# MCell: Molecules: setup information
	##########

	# Molecule names
	mol_names = ["s01","s11","s21","s31","s00","s10","s20","s30"]

	# Make all the molecule objects
	mol_obj_dict = {}
	# Go through all regions
	for sc_id, nsg in nsg_sc_dict.items():
		for isg in range(1,nsg+1):
			for mol in mol_names:
				mol_name = mol + "_sc_%02d_%02d_sg_%02d"%(sc_id[0],sc_id[1],isg)
				# print("Making molecule: " + str(mol_name))
				tmp_obj = MCell_Molecule(name = mol_name, type = mol, reg_id = (sc_id[0],sc_id[1],isg), dc = 0.0, surf_flag = True)
				mol_obj_dict[mol_name] = tmp_obj

	# Total molecules to release across all sections
	num_rel_dict_tot = {}
	num_rel_dict_tot[mol_names[0]] = 496
	num_rel_dict_tot[mol_names[1]] = 83
	num_rel_dict_tot[mol_names[2]] = 5
	num_rel_dict_tot[mol_names[3]] = 0
	num_rel_dict_tot[mol_names[4]] = 354
	num_rel_dict_tot[mol_names[5]] = 58
	num_rel_dict_tot[mol_names[6]] = 3
	num_rel_dict_tot[mol_names[7]] = 1

	# Factor
	factor = 10
	for key,val in num_rel_dict_tot.items():
		num_rel_dict_tot[key] = factor*val

	# Compute the number to release in each segment
	n_channel_in_reg_dict = calculate_geometry_info.compute_nrel(mol_obj_dict, num_rel_dict_tot, surf_area_dict)

	##########
	# MCell: Molecules: define in MCell
	##########

	# Go through all molecules
	for mol in mol_obj_dict.values():
		# Create the mols: world, name, DC, surface mol flag
		tmp_particle = m.create_species(world_main, mol.name, mol.dc, mol.surf_flag)
		mol.mcell_obj = tmp_particle

	##########
	# MCell: Reactions: setup information
	##########

	# Rate constants
	'''
	alpham = 143.434645141
	betam = 5665.90470964
	alphah = 95.7604961147
	betah = 25.9141573172
	'''

	# Resting potential (initial)
	v0 = -69.0

	# List of true rxns
	# Note: the reaction rates here "lambda" functions!
	rxn_list = [("s01", ";", "<-->", "s11", ";",lambda x: 3*f_alpham(x),lambda x: 1*f_betam(x)),
				("s11", ";", "<-->", "s21", ";",lambda x: 2*f_alpham(x),lambda x: 2*f_betam(x)),
				("s21", ";", "<-->", "s31", ";",lambda x: 1*f_alpham(x),lambda x: 3*f_betam(x)),
				("s00", ";", "<-->", "s10", ";",lambda x: 3*f_alpham(x),lambda x: 1*f_betam(x)),
				("s10", ";", "<-->", "s20", ";",lambda x: 2*f_alpham(x),lambda x: 2*f_betam(x)),
				("s20", ";", "<-->", "s30", ";",lambda x: 1*f_alpham(x),lambda x: 3*f_betam(x)),
				("s01", ";", "<-->", "s00", ";",lambda x: f_betah(x),lambda x: f_alphah(x)),
				("s11", ";", "<-->", "s10", ";",lambda x: f_betah(x),lambda x: f_alphah(x)),
				("s21", ";", "<-->", "s20", ";",lambda x: f_betah(x),lambda x: f_alphah(x)),
				("s31", ";", "<-->", "s30", ";",lambda x: f_betah(x),lambda x: f_alphah(x))]

	# List of molecule names
	mol_name_1_list = [item[0] for item in rxn_list]
	mol_name_2_list = [item[3] for item in rxn_list]

	##########
	# MCell: Create scene
	##########

	# Create a scene (?)
	scene_name = "Scene"
	scene = m.create_instance_object(world_main, scene_name)

	##########
	# MCell: Create the objects
	##########

	# Make the MCell obj
	mn_surface_obj = m.create_polygon_object(world_main, mn_surface_vert_list, mn_surface_el_list, scene, "MN_Surface")
	mn_segment_obj = m.create_polygon_object(world_main, mn_segment_vert_list, mn_segment_el_list, scene, "MN_Segment")
	
	##########
	# MCell: Create the surface regions AND surface classes
	##########

	# Create the surface regions AND surface classes
	mn_surface_surf_reg_dict = {}
	mn_segment_surf_reg_dict = {}
	mn_surface_surf_class_dict = {}
	mn_segment_surf_class_dict = {}
	
	# Surface
	for reg_name, reg_face_list in mn_surface_reg_dict.items():
		
		if not len(reg_face_list) > 0:
			print("ERROR HERE")
			print(reg_name)

		# Create the surface region
		mn_surface_surf_reg_dict[reg_name] = m.create_surface_region(world_main, mn_surface_obj, reg_face_list, reg_name)

		# Create the surface class
		surf_class_name = "SurfClass_" + reg_name
		mn_surface_surf_class_dict[surf_class_name] = m.create_surf_class(world_main, surf_class_name)
		m.mcell_assign_surf_class_to_region(mn_surface_surf_class_dict[surf_class_name], mn_surface_surf_reg_dict[reg_name])

	# Segment
	for reg_name, reg_face_list in mn_segment_reg_dict.items():

		# Create the surface region
		mn_segment_surf_reg_dict[reg_name] = m.create_surface_region(world_main, mn_segment_obj, reg_face_list, reg_name)

		# Create the surface class
		surf_class_name = "SurfClass_" + reg_name
		mn_segment_surf_class_dict[surf_class_name] = m.create_surf_class(world_main, surf_class_name)
		m.mcell_assign_surf_class_to_region(mn_segment_surf_class_dict[surf_class_name], mn_segment_surf_reg_dict[reg_name])

	##########
	# MCell: Reactions: define in MCell
	# Note that the reactions must be done AFTER defining the surface classes, since we have reactions that occur with surface classes
	##########

	# Dictionary of all real MCell reactions
	mcell_rxn_real_dict = {}
	# List of all inter-compartment reactions
	mcell_rxn_comp_list = []

	# Generate the reactions
	generate_rxn_list(mcell_rxn_real_dict, mcell_rxn_comp_list, mol_obj_dict, rxn_list, mol_names, v0, nsg_sc_dict, mn_segment_reg_dict, mn_segment_surf_class_dict)

	# Add reactions to MCell object
	# NOTE HERE
	# THERE IS SOME PROBLEM WITH REVERSIBLE SURFACE REACTIONS HERE
	# HENCE, DEFINE TWO REACTIONS....
	for rxns in mcell_rxn_real_dict.values():
		for rxn in rxns:
			m.create_reaction(world_main, rxn.reactant_list, rxn.product_list, rxn.rc_f, name=rxn.name+"_f")
			m.create_reaction(world_main, rxn.product_list, rxn.reactant_list, rxn.rc_b, name=rxn.name+"_b")

	for rxn in mcell_rxn_comp_list:
		m.create_reaction(world_main, rxn.reactant_list, rxn.product_list, rxn.rc_f, surf_class=rxn.surf_list, name=rxn.name+"_f")
		m.create_reaction(world_main, rxn.product_list, rxn.reactant_list, rxn.rc_b, surf_class=rxn.surf_list, name=rxn.name+"_b")

	##########
	# MCell: Release sites
	##########

	# Go through all molecules
	for mol in mol_obj_dict.values():
		rel_list = m.mcell_add_to_species_list(mol.mcell_obj, mol.surf_flag, 1, None)
		# Surface class name
		sc_name = "SurfClass_" + "sc_%02d_%02d_sg_%02d" % mol.reg_id + "_B_surf"
		# Create release obj
		throw_away = m.mcell_add_mol_release_to_surf_class(world_main, mn_surface_surf_class_dict[sc_name], rel_list, mol.nrel, 0, None)
	
	##########
	# MCell: Set up count statements
	##########

	# In order to count the number of s31, we need to have the following count statements
	print("> Setting up count statements")

	# Throwing the output away causes a memory error - the output of create_count MUST be assigned to a permanent object
	store_count_dict = {}

	for mol in mol_obj_dict.values():
		if mol.type == "s31":
			# Region name
			reg_name = "sc_%02d_%02d_sg_%02d" % mol.reg_id + "_B_surf"
			reg_sym = m.mcell_get_reg_sym(mn_surface_surf_reg_dict[reg_name])

			# create_count args: the world, the symbol to the object to count in, the mcell molecule, the filepath, the timestep at which to output
			store_count_dict[mol.name] = m.create_count(world_main, m.mcell_get_reg_sym(mn_surface_surf_reg_dict[reg_name]), mol.mcell_obj, "react_data/"+str(mol.name)+"_MN_Surface.dat", dt_M)

	print("> Set up count statements succesfully")

	##########
	# NEURON: Stimulation
	##########

	#'''
	fsec = nrn_sc_dict[(1,2)]
	stim = h.IClamp(fsec(0.0))
	stim.delay = 0.0 # [ms] delay
	stim.dur = 0.1  # [ms] duration
	stim.amp = 0.1   # [nA] amplitude
	#'''

	##########
	# MCell: Initialize the MCell simulation
	##########

	print("> Initializing MCell simulation")

	m.mcell_init_simulation(world_main)
	m.mcell_init_output(world_main)

	print("> Finished init MCell simulation")

	##########
	# NEURON: Initialize membrane potential
	##########

	print("> Initializing NEURON membrane potential")

	# Initial voltage value
	h.finitialize(v0)

	print("> Finished init NEURON membrane potential")

	##########
	# Run
	##########

	print("> Starting Simulation")

	# Go through all NEURON iterations
	for i_iter_N in range(0,n_iter_N):

		# print("> NEURON iteration: " + str(i_iter_N) + " / " + str(n_iter_N))

		##########
		##########
		# MCell
		##########
		##########

		# First set the rates as the correct values in all segments at these voltage values

		# Dictionary of the fraction of open channels in each region at half a timestep
		frac_open_dict = {}

		###
		# Set the rate constants based on the current voltage values in each segment
		###

		# Go through all real reactions
		for rxns in mcell_rxn_real_dict.values():
			for rxn in rxns:

				# Voltage in the region this reaction occurs in
				v_reg = nrn_sg_dict[rxn.reg_id].v

				# Compute the new RC
				# Forward rxn
				if rxn.f_rc_f != None:
					rxn.rc_f = rxn.f_rc_f(v_reg)
					m.mcell_modify_rate_constant(world_main, rxn.name + "_f", rxn.rc_f)
				# Backward
				if rxn.f_rc_b != None:
					rxn.rc_b = rxn.f_rc_b(v_reg)
					m.mcell_modify_rate_constant(world_main, rxn.name + "_b", rxn.rc_b)

		# print("> > Finished modifying rate constants")

		###
		# Run MCell
		###

		# Run the "normal" world (normal voltages) for full a timestep
		output_freq = 1 # Output at every timestep, rather than at every "output_freq"-th
		for i_mstep in range(0,n_steps_M):
			# args: world, output frequency, verbose: 0 = normal output; 1 = quiet
			m.mcell_run_iteration(world_main, output_freq, 1)

			# Count the number of open channels in each region
			for sec_id,nseg in nsg_sc_dict.items():
				for i_seg in range(1,nseg+1):

					# Region id and name
					reg_id = (sec_id[0],sec_id[1],i_seg)
					reg_name = "sc_%02d_%02d_sg_%02d" % reg_id

					# Count the number of open channels
					region_name = reg_name + "_B_surf"
					open_channel_count = m.mcell_get_count("s31_"+reg_name, "%s.%s,%s" % (scene_name, "MN_Surface", region_name), world_main)

					if reg_id in frac_open_dict:
						frac_open_dict[reg_id] += open_channel_count
					else:
						frac_open_dict[reg_id] = open_channel_count

		# Convert open channels to a fraction
		for key,val in frac_open_dict.items():

			# Store the current in this region
			frac_open_dict[key] = (1.0 * val) / n_channel_in_reg_dict[key] / n_steps_M

		# print("frac open: " + str(list(frac_open_dict.values())[0]))

		# print("> > MCell: Ran full timestep")

		##########
		##########
		# NEURON
		##########
		##########

		###
		# Set the fraction of open channels in neuron
		###

		# print("> > Fraction of open channels:")
		# print(frac_open_dict)

		# Go through all the segments
		for reg_id, sg in nrn_sg_dict.items():
			for mech in sg:
				if mech.name() == "MCell_na":
					# Set the currents
					mech.frac_open = frac_open_dict[reg_id]
		
		###
		# Advance NEURON by one timestep
		###

		# Advance neuron
		h.fadvance()

		##########
		##########
		# Print voltages
		##########
		##########

		print("> > Finished; time: " + str(dt_N * i_iter_N) + " ; voltage in: " + str((1,2,1)) + " = " + str(nrn_sg_dict[(1,2,1)].v))

		##########
		##########
		# Write voltage data to file
		##########
		##########

		f = open(DIR+"v_data/v_%03d.txt"%i_iter_N,'w')
		FIRST_FLAG = True
		for reg_id, sg in nrn_sg_dict.items():
			if not FIRST_FLAG:
				f.write("\n")
			else:
				FIRST_FLAG = False
			f.write(str(reg_id[0]) + " " + str(reg_id[1]) + " " + str(reg_id[2]) + " " + str(sg.v))
		f.close()


	m.mcell_flush_data(world_main)
	m.mcell_print_final_warnings(world_main)
	m.mcell_print_final_statistics(world_main)

