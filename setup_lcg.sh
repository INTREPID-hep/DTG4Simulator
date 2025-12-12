#!/bin/bash

# --- Configuration ---
# Preferred LCG Version. 
# Note: Newer LCG versions (like 108) might not exist for older OSs (like CentOS 7).
LCG_VERSION_EL9="LCG_108"   # Contains Geant4 11.3.2
LCG_VERSION_C7="LCG_105"    # Contains Geant4 11.2.0 (Last one supporting CentOS 7)

# --- Automatic Detection ---
# Only print messages if we are in an interactive shell or explicitly asked
if [[ $- == *i* ]]; then
    VERBOSE=1
else
    VERBOSE=0
fi

[ $VERBOSE -eq 1 ] && echo "Detecting system..."
ARCH=$(uname -m)
OS_TAG="unknown"

if [ -f /etc/os-release ]; then
    source /etc/os-release
    if [[ "$PLATFORM_ID" == "platform:el9" ]]; then
        OS_TAG="el9"
        TARGET_LCG=$LCG_VERSION_EL9
    elif [[ "$ID" == "centos" && "$VERSION_ID" == "7" ]]; then
        OS_TAG="centos7"
        TARGET_LCG=$LCG_VERSION_C7
    fi
fi

if [ "$OS_TAG" == "unknown" ]; then
    echo "Error: Could not detect a supported OS (EL9 or CentOS 7)."
    return 1 2>/dev/null || exit 1
fi

[ $VERBOSE -eq 1 ] && echo "Detected OS: $OS_TAG"
[ $VERBOSE -eq 1 ] && echo "Target LCG Version: $TARGET_LCG"

# Find the best compiler view (prefer gcc, optimized)
# We search for a directory matching the pattern in the specific LCG version
VIEW_DIR=$(find /cvmfs/sft.cern.ch/lcg/views/$TARGET_LCG -maxdepth 1 -name "${ARCH}-${OS_TAG}-gcc*-opt" | sort -V | tail -n 1)

if [ -z "$VIEW_DIR" ]; then
    echo "Error: Could not find a suitable LCG view in /cvmfs/sft.cern.ch/lcg/views/$TARGET_LCG"
    return 1 2>/dev/null || exit 1
fi

SETUP_SCRIPT="$VIEW_DIR/setup.sh"

[ $VERBOSE -eq 1 ] && echo "Loading LCG View: $VIEW_DIR"
source "$SETUP_SCRIPT"
