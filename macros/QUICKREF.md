# DTSim Macro Quick Reference

## Running Simulations

### Interactive Mode (Visualization)
```bash
# Default - loads macros/interactive/vis.mac
./exampleDTSim

# With magnetic field
./exampleDTSim -B
```

### Batch Mode
```bash
# Basic batch run (10 events)
./exampleDTSim -m macros/batch/run.mac

# High-pT muons (1 TeV, 100 events, 4 threads)
./exampleDTSim -m macros/batch/run_highpt.mac

# Low-pT mixed particles (10 GeV randomized, 100 events)
./exampleDTSim -m macros/batch/run_lowpt.mac
```

## Customizing Settings

### Change Default Particle
Edit `macros/settings/primary.mac`:
```plaintext
/gun/particle pi+              # Change particle type
/DTSim/generator/momentum 500 GeV  # Change momentum
```

### Change Visualization Colors
Edit `macros/visualization/vis_trajectories.mac`:
```plaintext
/vis/modeling/trajectories/drawByParticleID-0/set mu+ red  # Change muon color
```

### Change Camera View
Edit `macros/visualization/vis_common.mac`:
```plaintext
/vis/viewer/set/viewpointThetaPhi 90. 0.  # Side view
/vis/viewer/zoom 2.0                       # Zoom in
```

## Creating New Macros

### Template for Batch Run
```plaintext
# my_custom_run.mac
/control/execute macros/settings/verbosity.mac
/control/execute macros/settings/physics.mac
/run/numberOfThreads 4
/run/initialize

# Custom settings
/DTSim/generator/momentum 100 GeV
/gun/particle proton

/run/beamOn 1000
```

### Template for Interactive Run with Custom Visualization
```plaintext
# my_custom_vis.mac
/control/execute macros/settings/verbosity.mac
/control/execute macros/settings/physics.mac
/run/initialize

# Custom visualization (instead of using vis_common.mac)
/vis/open OGL 800x800-0+0
/vis/viewer/set/viewpointThetaPhi 45. 45.
/vis/drawVolume

# Still use standard components
/control/execute macros/visualization/vis_scene.mac
/control/execute macros/visualization/vis_trajectories.mac
/control/execute macros/settings/primary.mac

/vis/viewer/set/autoRefresh true
/vis/viewer/flush
```

## Common Commands

### Particle Gun
```plaintext
/gun/particle <particle>       # e.g., mu+, pi+, proton, gamma
/gun/energy 10 GeV
/gun/position 0 0 -0.5 m
```

### DTSim Generator
```plaintext
/DTSim/generator/momentum 300 GeV
/DTSim/generator/sigmaMomentum 50 MeV
/DTSim/generator/sigmaAngle 2 deg
/DTSim/generator/randomizePrimary true
```

### Visualization
```plaintext
/vis/viewer/set/style surface         # or wireframe
/vis/viewer/zoom 1.5
/vis/viewer/refresh
/vis/scene/endOfEventAction accumulate  # or refresh
```

### Run Control
```plaintext
/run/numberOfThreads 8
/run/beamOn 100
/run/verbose 2
```

## File Locations

- **Reusable Settings**: `macros/settings/`
- **Visualization Components**: `macros/visualization/`
- **Interactive Macros**: `macros/interactive/`
- **Batch Macros**: `macros/batch/`
- **Original Macros (backup)**: `macros/old/`
