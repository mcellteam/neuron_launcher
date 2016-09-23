# Neuron Launcher
Set of tools in Blender to create models for MCell/NEURON hybrid simulations

## Requirements

Blender, CellBlender, and possible MeshPy based on the usage scenario.

## Installation

The default method for installing Blender add-ons suffices:

	Download/Clone the repo as a compressed zip file
	In Blender, under User Preferences > Add-ons click "Install from File"
	Remember to activate the Plug-In

## Model Creation

The tools included in this add-on are designed to create models for hybrid MCell/NEURON simulations, starting from a cable model in the form of an SWC file. Included are Blender scripts for modifying the meshes and cable model, as well as the software SWC Mesher, which extrapolates cable models (SWC files) to create a surface mesh.

* **[Description](readme_files/description)**
* **[SWC Mesher](https://github.com/mcellteam/swc_mesher)**
* **[Closing an open mesh](readme_files/closing_mesh)**
* **[Assigning surface regions](readme_files/assigning_surface_regions)**
* **[Create volumetric compartments from the surface mesh](readme_files/creating_compartments)**
* **[Visualising exploding sections](readme_files/visualising_exploding_sections)**
* **[Source file list](readme_files/source_file_list)**

![Example](readme_files/figures/example.jpg?raw=true "Example")