import math

import random

# Function to get the number of segments per section
def compute_nsg_sc(reg_dict):

    nsg_sc_dict = {}

    for reg in reg_dict.keys():
        if len(reg) == 21 and reg[-6:] == "B_surf":
            sc_id = (int(reg[3:5]),int(reg[6:8]))
            sg_id = int(reg[12:14])
            if not sc_id in nsg_sc_dict:
                nsg_sc_dict[sc_id] = sg_id
            else:
                if sg_id > nsg_sc_dict[sc_id]:
                    nsg_sc_dict[sc_id] = sg_id

    return nsg_sc_dict

# Function to compute surface areas for each regions
def compute_surface_areas(vert_list, el_list, reg_dict):

    surf_area_dict = {}

    # Go through all regions
    for reg_name, reg_f_list in reg_dict.items():

        sa = 0.0

        # Go through all faces
        for f in reg_f_list:

            vs = el_list[f]
            cos = [list(vert_list[item]) for item in vs]
            d1 = [cos[1][i] - cos[0][i] for i in [0,1,2]]
            d2 = [cos[2][i] - cos[0][i] for i in [0,1,2]]
            crss = [d1[1] * d2[2] - d1[2] * d2[1],
                    d1[2] * d2[0] - d1[0] * d2[2],
                    d1[0] * d2[1] - d1[1] * d2[0]]
            sa += 0.5*((crss[0]**2 + crss[1]**2 + crss[2]**2)**(1/2))

        # Convert the region name into ids
        reg_ids = (int(reg_name[3:5]),int(reg_name[6:8]),int(reg_name[12:14]))

        surf_area_dict[reg_ids] = sa

    return surf_area_dict

# Function to compute the number of molecules to release into each region
# Input: num_rel_dict should be a dictionary of molecule name => number to release (in all regions total)
def compute_nrel(mol_obj_dict, num_rel_dict, surf_area_dict):

    # Convert a surface area dictionary to a fractional area dictionary
    surf_area_tot = sum(list(surf_area_dict.values()))
    frac_area_dict = {}
    for reg, val in surf_area_dict.items():
        frac_area_dict[reg] = val/surf_area_tot

    # Dictionary of region name to molecule name to number to release
    rel_dict = {}

    # Go through all regions and determine how many of each molecule to release
    nremaining_dict = num_rel_dict
    for reg, frac in frac_area_dict.items():
        rel_dict[reg] = {}

        # Go through every molecule
        for mol, nrel_tot in num_rel_dict.items():

            nrel = int(frac*nrel_tot)
            rel_dict[reg][mol] = nrel
            nremaining_dict[mol] -= nrel

    # There are some missing; assign these randomly
    random.seed(2)
    reg_names = list(rel_dict.keys())
    for mol, nrel in nremaining_dict.items():
        for i_rel in range(0,nrel):
            i_rand = random.randint(0,len(reg_names)-1)
            rel_dict[reg_names[i_rand]][mol] += 1

    # Assign the number to release to mcell objects
    for mol in mol_obj_dict.values():
        mol.nrel = rel_dict[mol.reg_id][mol.type]

    # Also return the total number of channels per region
    n_channel_in_reg_dict = {}
    for reg_id, rel_mol_dict in rel_dict.items():
        nrel_reg = 0
        for nrel_mol in rel_mol_dict.values():
            nrel_reg += nrel_mol

        n_channel_in_reg_dict[reg_id] = nrel_reg

    return n_channel_in_reg_dict


