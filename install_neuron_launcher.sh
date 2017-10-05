#!/usr/bin/env bash

# neuron_launcher will be installed here. This should typically be a link to the desired location on your platform.
# For example, for a link pointing into a MacOSX bundle: ln -s /Applications/Blender-2.78c-CellBlender/blender.app/Contents/Resources/2.78/scripts/addons/ ~/my_blender_addons_link

INSTALL_DIR=~/my_blender_addons_link/

# Install neuron launcher
if [ -d "$INSTALL_DIR/neuron_launcher" ]; then
	rm -r "$INSTALL_DIR/neuron_launcher"
fi
mkdir "$INSTALL_DIR/neuron_launcher"
cp src/* "$INSTALL_DIR/neuron_launcher/"
