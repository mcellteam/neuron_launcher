
INSTALL_DIR = ~/my_blender_addons_link

SHELL = /bin/sh

SOURCES = ./neuron_launcher/check_connectivity.py ./neuron_launcher/check_double_assignments.py ./neuron_launcher/check_intersecting_faces.py ./neuron_launcher/check_overlapping_faces.py ./neuron_launcher/check_unassigned_faces.py ./neuron_launcher/close_open_caps.py ./neuron_launcher/color_regions.py ./neuron_launcher/compartmentize_cyl.py ./neuron_launcher/compartmentize_sc_only.py ./neuron_launcher/compartmentize_tet.py ./neuron_launcher/explode.py ./neuron_launcher/__init__.py ./neuron_launcher/neuron_launcher_gui.py ./neuron_launcher/regions_to_compartments.py ./neuron_launcher/select_faces.py ./neuron_launcher/surface_sections.py ./neuron_launcher/swc_mesher.py ./neuron_launcher/timeline_voltage.py

ZIPFILES = $(SOURCES)

ZIPOPTS = -X -0 -D -o


all: neuron_launcher neuron_launcher.zip


neuron_launcher:
	ln -s . neuron_launcher


neuron_launcher.zip: $(SOURCES)
	@echo Updating neuron_launcher.zip
	@echo Sources = $(SOURCES)
	@zip $(ZIPOPTS) neuron_launcher.zip $(ZIPFILES)


clean:
	rm -f neuron_launcher.zip


install: neuron_launcher.zip
	@ mkdir -p $(INSTALL_DIR)
	@ unzip -o neuron_launcher.zip -d $(INSTALL_DIR); \

