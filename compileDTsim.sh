#!/bin/bash

# Source the environment setup script
# This ensures we have the correct compiler and Geant4 version
echo $(dirname "$0")/setup_lcg.sh

# --- Build Steps ---
echo "Creating build directory"
mkdir -p DTSim_build
cd DTSim_build

echo "Configuring..."
cmake ../

echo "Compiling..."
make -j4

