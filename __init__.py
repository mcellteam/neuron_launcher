# -*- coding: utf-8 -*-
"""
Created on Tue Sept 20 00:00:00 2016

@author: Oliver K. Ernst <oernst@salk.edu>
"""

bl_info = {
    "name": "Neuron Launcher",
    "description": "Tools for MCell-NEURON hybrid model construction",
    "author": "Oliver K. Ernst, Tom Bartol, Bob Kuczewski",
    "version": (0,1,0),
    "blender": (2, 7, 7),
    "location": "View3D > Add > Mesh",
    "warning": "",
    "wiki_url" : "http://salk.edu",
    "license": "GPL v2",
    "category": "Mesh"}

# -------------------------------------------------------------------------- 
# ***** BEGIN GPL LICENSE BLOCK ***** 
# 
# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software Foundation, 
# Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
# 
# ***** END GPL LICENCE BLOCK ***** 
# -------------------------------------------------------------------------- 

if "bpy" in locals():
    print("Reloading Neuron Launcher")
    import imp
    imp.reload(swc_mesher)
    imp.reload(neuron_launcher_gui)
    imp.reload(close_open_caps)
    imp.reload(surface_sections)
    imp.reload(check_bordering_surfaces)
    imp.reload(check_double_assignments)
    imp.reload(color_regions)
    imp.reload(compartmentize_tet)
    imp.reload(compartmentize_cyl)
    imp.reload(compartmentize_sc_only)
    imp.reload(explode)
    imp.reload(regions_to_compartments)

else:
    print("Importing Neuron Launcher")
    from . import swc_mesher
    from . import neuron_launcher_gui
    from . import close_open_caps
    from . import surface_sections
    from . import check_bordering_surfaces
    from . import check_double_assignments
    from . import color_regions
    from . import compartmentize_tet
    from . import compartmentize_cyl
    from . import compartmentize_sc_only
    from . import explode
    from . import regions_to_compartments

# General import
import bpy
import sys
import os

def register():
    bpy.utils.register_module(__name__)
    
    bpy.types.Scene.make_neuron_meta = bpy.props.PointerProperty(type=swc_mesher.MakeNeuronMetaPropGroup)
    bpy.types.Scene.nrnlauncher = bpy.props.PointerProperty(type=neuron_launcher_gui.NeuronLauncherPropGroup)

    print("Neuron Launcher registered")



def unregister():
    
    bpy.utils.unregister_module(__name__)

    print("Neuron Launcher unregistered")

# ?
if __name__ == '__main__':
    register()
