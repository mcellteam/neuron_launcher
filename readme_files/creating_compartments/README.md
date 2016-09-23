# Creating compartments

It is not sufficient to assign surface patches to different sections/segments in the cable model; volumetric compartments must be assigned as well. This means that dividing planes must be created through the volume that separate different compartments in MCell, and then each compartment must be associated with a different segment in the NEURON cable model.

Two methods have been developed to deal with this, described below. Future work will also allow combinations of these methods.

## Tetrahedralization method

The most robust method is based on tetrahedralization of the volume. The entire volume is tetrahedralized using MeshPy, a Python interface to TetGen. Afterwards, individual tets are assigned to different segments in the cable model forming separate compartments. Finally, dividing planes between the compartments are created. This method is powerful, but very slow, both due to the tetrahedralization and the computation involved in assigning tets to segments. Optimally, it should be used for non-trivial geometries, e.g. for non-cable-like structures such as the soma.

The figure below shows the starting point for this function - MCell regions must have been assigned as described [here](../assigning_surface_regions).

![Starting Point](../figures/creating_compartments_1.jpg?raw=true "Starting Point")

It is best to monitor the progress in the terminal, as it may be slow.

The output generates two new objects: a surface object which contains surface MCell regions for every segment, and a segment object which contains the dividing planes between the compartments, shown below.

![Output](../figures/creating_compartments_2.jpg?raw=true "Output")
![Dividing Planes](../figures/creating_compartments_3.jpg?raw=true "Dividing Planes")


## Cylinder method

A faster method makes use of the cylindrical-like structure of the surface mesh to create logical subdivisions. Currently, this is only implemented for dividing sections, not further subdivisions into segments. Note that this will work well for tube-like structures, but not as well for non-trivial geometries such as the soma.

The starting point is shown below.

![Starting Point Sections Only](../figures/creating_compartments_4.jpg?raw=true "Starting Point Sections Only")

Pressing the "Compartmentize Sections Only" button gives the output shown below.

![Output Sections Only](../figures/creating_compartments_5.jpg?raw=true "Output Sections Only")
![Dividing Planes Sections Only](../figures/creating_compartments_6.jpg?raw=true "Dividing Planes Sections Only")

# Converting MCell regions into separate objects

A further tool is included to convert the output of the compartmentalization from MCell regions into separate objects. This is particularly useful for visualisation.

Start by selecting both of the output objects from either of the two compartmentize methods described above, i.e. the "_Segment" object which contains the dividing planes between compartments, and the "_Surface" object which contains the surface regions.

![Prepare to make compartments](../figures/creating_compartments_6.jpg?raw=true "Prepare to make compartments")

The "Convert Regions to Compartments" will generate separate objects for each compartment. These may either be segment compartments, or just the larger section compartments, depending on what is defined in the surface and segment object regions.

The output is shown below, with several of the compartments hidden for illustration, and one selected.

![Output Compartments](../figures/creating_compartments_7.jpg?raw=true "Output Compartments")



