#!/usr/bin/env bash

# Find the system type and Python version
if [ "$(uname)" == "Darwin" ]; then
    # Mac OS X
    INSTALL_DIR=~/Library/Application\ Support/Blender/2.78/scripts/addons
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # GNU/Linux
    INSTALL_DIR=~/.config/blender/2.78/scripts/addons
else
	echo "Error: cannot determine system - check install script."
	return
fi

# Install neuron launcher
if [ -d "$INSTALL_DIR/neuron_launcher" ]; then
	rm -r "$INSTALL_DIR/neuron_launcher"
fi
mkdir "$INSTALL_DIR/neuron_launcher"
cp src/* "$INSTALL_DIR/neuron_launcher/"