#!/bin/bash

# Source the environment setup script
# This ensures we have the correct compiler and Geant4 version
source $(dirname "$0")/setup_lcg.sh

# --- Build Steps ---
echo "Creating build directory"
mkdir -p build
cd build

echo "Configuring..."
cmake ../

echo "Compiling..."
make -j4

