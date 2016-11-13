# Neuron Launcher
Set of tools in Blender to create models for MCell/NEURON hybrid simulations

## Requirements

Blender, CellBlender, and MeshPy based on the usage scenario. MeshPy is included in this repo - CellBlender can be found [here](https://github.com/mcellteam/cellblender).

## Installation

NeuronLauncher requires [MeshPy](https://mathema.tician.de/software/meshpy/), which is contained here for convenience along with a simple install script.

	Download/Clone the repo
	Install MeshPy if necessary by running:
	sh install_meshpy.sh
	Install NeuronLauncher by running:
	sh install_neuron_launcher.sh

## Model Creation

The tools included in this add-on are designed to create models for hybrid MCell/NEURON simulations, starting from a cable model in the form of an SWC file. Included are Blender scripts for modifying the meshes and cable model, as well as the software SWC Mesher, which extrapolates cable models (SWC files) to create a surface mesh.

* **[Description](readme_files/description)**
* **[SWC Mesher](https://github.com/mcellteam/swc_mesher)**
* **[Closing an open mesh](readme_files/closing_mesh)**
* **[Assigning surface regions](readme_files/assigning_surface_regions)**
* **[Create volumetric compartments from the surface mesh](readme_files/creating_compartments)**
* **[Visualising exploding sections](readme_files/exploding_sections)**
* **[Visualising Materials and Colors](readme_files/visualising)**
* **[Source file list](readme_files/source_file_list)**

![Example](readme_files/figures/example.jpg?raw=true "Example")
