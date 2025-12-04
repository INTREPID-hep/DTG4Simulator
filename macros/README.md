# DTSim Macros Organization

This directory contains organized macro files for running DTG4Simulator simulations.

## Directory Structure

```
macros/
├── settings/          # Reusable configuration blocks
├── visualization/     # Modular visualization components
├── interactive/       # Interactive mode macros
├── batch/            # Batch/production run macros
└── old/              # Original macros (backup)
```

## Settings Macros (`settings/`)

Reusable configuration blocks that can be included in any macro:

- **`verbosity.mac`** - Standard verbose levels for all managers
- **`physics.mac`** - Physics list settings, cuts, threading
- **`primary.mac`** - Default primary particle configuration

## Visualization Macros (`visualization/`)

Modular visualization components for building custom views:

- **`vis_common.mac`** - Opens viewer, sets viewport, style, zoom
- **`vis_trajectories.mac`** - Trajectory modeling, filters, colors by particle ID
- **`vis_scene.mac`** - Scene elements (axes, hits, magnetic field, event ID)

## Interactive Macros (`interactive/`)

User-facing macros for interactive visualization:

- **`vis.mac`** - Main interactive visualization macro (auto-loaded)
- **`gui.mac`** - GUI menu bar configuration
- **`icons.mac`** - Icon toolbar configuration

## Batch Macros (`batch/`)

Production/batch running macros:

- **`run.mac`** - Basic batch execution
- **`run_highpt.mac`** - High-pT muon run (1 TeV, 100 events)
- **`run_lowpt.mac`** - Low-pT mixed particle run (10 GeV, randomized)

## Usage Examples

### Interactive Mode (with visualization)
```bash
./exampleDTSim
# Automatically loads macros/interactive/vis.mac
```

### Batch Mode
```bash
./exampleDTSim -m macros/batch/run.mac
./exampleDTSim -m macros/batch/run_highpt.mac
./exampleDTSim -m macros/batch/run_lowpt.mac
```

### Creating Custom Macros

You can easily create new macros by including existing components:

```plaintext
# my_custom_run.mac
/control/execute macros/settings/verbosity.mac
/control/execute macros/settings/physics.mac
/run/initialize

# Custom primary settings
/DTSim/generator/momentum 500 GeV
/gun/particle pi+

/run/beamOn 50
```

### Modifying Visualization

To change visualization settings globally, edit the files in `visualization/`:
- Change particle colors → edit `vis_trajectories.mac`
- Change camera angle/zoom → edit `vis_common.mac`
- Add/remove scene elements → edit `vis_scene.mac`

All macros using these components will automatically pick up the changes.

## Benefits

✅ **Reusability** - Common settings shared across all macros  
✅ **Maintainability** - Change once, applies everywhere  
✅ **Clarity** - Each file has a single, clear purpose  
✅ **Flexibility** - Easy to create new scenarios by mixing components
