#!/usr/bin/env bash

# Find the system type and Python version
if [ "$(uname)" == "Darwin" ]; then
    # Mac OS X
    PY_DIR=/Applications/blender.app/Contents/Resources/2.78/python/bin/python3.5m
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # GNU/Linux
    PY_DIR=~/bin/blender/2.78/python/bin/python3.5m
else
	echo "Error: cannot determine system - check install script."
	return
fi

# meshpy configure
cd MeshPy-2016.1.2
$PY_DIR configure.py
$PY_DIR setup.py install
