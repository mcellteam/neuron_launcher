# Assigning surface regions

Surface regions are defined in CellBlender as MCell regions. To create a correspondence between the surface mesh and the cable model, each face must be assigned to a section in the cable model (to a straight line segment). This may be done using the "Make MCell surface regions" tool.

First, the SWC file for the object must be defined. Select the object, and under "Mesh objects and corresponding SWC files", add the object using the plus icon. Then click "Set SWC file" to set the SWC file corresponding to the object, as shown below.

![SWC assigned](../figures/assigning_surface_regions_1.jpg?raw=true "SWC assigned")

Next, the surface regions are calculated using the "Make MCell surface regions tool". Note that some faces may not be assigned if a "good" correspondence to a cable section cannot be found. Monitor the output of Blender in terminal to ensure that all faces are assigned; remaining faces must be assigned by hand.

The output of the function is that:
(1) The object is added to CellBlender, if not already
(2) MCell regions are created for each section. The naming convention is "sc_##_##", where "sc" stands for section, and the two numbers are the points in the SWC file which this section connects, with the lower number first and the index starting at 1.

![Surface Regions](../figures/assigning_surface_regions_2.jpg?raw=true "Surface Regions")

Many problems may occur as a result of poor assignments. In particular, the most devastating is if two faces border each-other from different sections, but these sections are not neighbors in the SWC file (i.e. in the cable model). This will be documented further later.

Two functions currently exist to monitor the quality of the assignments:
* "Check bordering surface assignments" checks the above problem. Note: SLOW!
* "Check for double assigned faces" checks that every face is assigned to only one MCell region.
