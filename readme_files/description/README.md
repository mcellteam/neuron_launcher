# Description

The purpose of these tools is for hybrid MCell/NEURON model creation.

* MCell models consist of a surface in 3D space, in which reaction-diffusion systems occur, and
* NEURON models consist of cable models, consisting of straight line pieces called "sections", which are further subdivided in "segments" (note that sections are physical objects, while segments are computational constructs to define the resolution at which to compute the voltage along the cable).

These tools prepare hybrid models by establishing correspondences between these two model types.

**The current starting point is a cable model in the form of an SWC file** (e.g. from the NeuroMorpho database, or ModelDB, or your favorite repo.). Eventually, methods for starting from SEM data will be developed.

The workflow proceeds as follows:
* **Import the cable model to Blender:** this may be done using the included SWC Mesher tool.
* **Extrapolate a surface mesh from the cable model:** this is also done using the SWC Mesher tool.
* **Cut the ROI out of the surface mesh and cable model:** this must currently be done by hand in Blender.
* **Close the surface mesh:** this may be done using the ["close an open mesh"](../closing_mesh) tool.
* **Assign surface regions:** this may be done using the ["assign surface regions"](../assigning_surface_regions) tool
* **Create volumetric compartments from the surface mesh:** this may be done using the ["compartmentalize"](../creating_compartments) tools.
* **Visualisation:** some tools for visualisation are included under ["visualising exploding sections"](../visualising).