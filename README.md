# DTG4Simulator

Simulación Geant4 de cámaras Drift Tubes (DT) del CMS para estudios de pattern recognition con muones y showers electromagnéticos.

## Requisitos

- **Geant4** 11.2+ con Qt/OpenGL
- **CMake** 3.16+
- **Python** 3.9+ con **mplDTs**: `pip install git+https://github.com/DanielEstrada971102/mplDTs.git@v2.2.0-beta`
- **ROOT** 6.x

## Instalación y Uso

```bash
# Clonar y compilar
git clone https://github.com/INTREPID-hep/DTG4Simulator.git
cd DTG4Simulator
mkdir build && cd build
cmake .. && make -j$(nproc)

# Ejecutar
./exampleDTSim                    # Modo interactivo con visualización
./exampleDTSim -b -m run.mac      # Modo batch
./exampleDTSim -B -m run.mac      # Con campo magnético (-2T en yoke)
```

## Salida de Datos

Archivos ROOT `DTG4Simulation_{runID}.root` con NTuple `DTG4SimNTuple/DTG4Tree` (14 columnas):
- **Identificación**: eventNumber, PDG, q, wheel, sector, station, superlayer, layer, cell
- **Posición**: xlocal, ylocal, zlocal [mm]
- **Física**: timewithdrift [ns], edep [MeV]

## Documentación Técnica

- **[Physics List](docs/physics.md)**: FTFP_BERT, campo magnético, control de tracking
- **[Geometría](docs/geometry.md)**: Generación desde Python con mplDTs, formato Text Geometry
- **[Primary Generator](docs/primary_generation.md)**: Configuración del particle gun, comandos UI
- **[Sensitive Detectors](docs/sensitive_detectors.md)**: DriftCellSD, decodificación de CellID, drift time
- **[Análisis](docs/analysis.md)**: Estructura del NTuple
## Modificar Geometría

```bash
cd geometry/
# Editar estaciones en generate_geometry_data.py
python generate_geometry_data.py
cd ../build && cmake ..  # Copiar nuevos .tg
```

## Producción Masiva

```bash
python submitJobs.py  # Genera y ejecuta múltiples jobs
hadd merged.root DTG4Simulation_*.root  # Combinar salidas
```

## Estructura del Proyecto

```
DTG4Simulator/
├── exampleDTSim.cc          # Main
├── include/src/             # Headers e implementaciones
├── geometry/                # Scripts Python + archivos .tg
├── macros/                  # vis.mac, run.mac, gui.mac
├── docs/                    # Documentación técnica detallada
└── build/                   # Compilación
```

## Licencia

Geant4 BSD License - Ver [LICENSE](LICENSE)
