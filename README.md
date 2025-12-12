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
./exampleDTSim                           # Modo interactivo con visualización
./exampleDTSim macros/batch/run.mac      # Modo batch con macro específico
```

## Salida de Datos

Archivos ROOT `DTG4Simulation_{runID}.root` con NTuple `DTG4SimNTuple/DTG4Tree`.
La estructura es **basada en vectores** (una fila por evento) y contiene cuatro colecciones principales:

- **Gen**: Información del generador (`gen_pt`, `gen_eta`, `gen_pdgId`...)
- **SimHits**: Hits de Geant4 (`simHit_wheel`, `simHit_xlocal`, `simHit_time`...)
- **Digis**: Señales digitalizadas (`digi_wheel`, `digi_TDC`...)
- **Muon Segments**: Trayectorias de muones por estación (`seg_wheel`, `seg_localPosX`, `seg_localDirX`...)

Para la lista completa de variables, ver [docs/analysis.md](docs/analysis.md).

## Configuración Runtime

Todos los parámetros del detector, física y digitalización son configurables mediante macros sin recompilar:

### Detector y Geometría (ANTES de /run/initialize)
```bash
/DTSim/detector/setGeometryFile geometry/geometry_concentrator.tg
/DTSim/detector/enableDriftSD true
/DTSim/detector/enableStationSD true
/DTSim/detector/useBField true
/DTSim/detector/BField/setGlobal 0 0 0 tesla
/DTSim/detector/BField/setYoke 0 0 -2.0 tesla
```

### Física del Detector (DESPUÉS de /run/initialize)
```bash
/DTSim/cellSD/setDriftVelocity 0.054 mm/ns    # 54 μm/ns
/DTSim/cellSD/setMinEnergy 26.6 eV
/DTSim/cellSD/setBarrierEnergy 2.1 keV
/DTSim/cellSD/setWallLoss 1.0 keV
/DTSim/cellSD/enableElectrostaticConfinement true
/DTSim/cellSD/enableWallCrossing true
```

### Digitalización (DESPUÉS de /run/initialize)
```bash
/DTSim/digitizer/setEfficiency 1.0
/DTSim/digitizer/setTimeResolution 2.0 ns
/DTSim/digitizer/setTDCResolution 0.78125 ns
```

Ver archivos en `macros/settings/` para ejemplos completos. Consultar **[docs/configuration.md](docs/configuration.md)** para referencia completa de todos los comandos disponibles.

## Documentación Técnica

- **[Configuración Runtime](docs/configuration.md)**: Referencia completa de comandos UI messenger
- **[Physics List](docs/physics.md)**: FTFP_BERT, campo magnético, control de tracking
- **[Geometría](docs/geometry.md)**: Generación desde Python con mplDTs, formato Text Geometry
- **[Primary Generator](docs/primary_generation.md)**: Configuración del particle gun, comandos UI
- **[Sensitive Detectors](docs/sensitive_detectors.md)**: DriftCellSD, decodificación de CellID, drift time
- **[Digitalización](docs/digitization.md)**: Simulación de electrónica, eficiencia, resolución temporal
- **[Análisis](docs/analysis.md)**: Estructura del NTuple
## Modificar Geometría

```bash
cd geometry/
# Editar estaciones en generate_geometry_data.py
python generate_geometry_data.py
cd ../build && cmake ..  # Copiar nuevos .tg
```

## Producción Masiva (HTCondor)

El script `submitJobs.py` automatiza la generación de macros y el envío de trabajos a HTCondor.

1.  **Configurar**: Editar `submitJobs.py` para definir los datasets, número de jobs, comandos específicos y directorio de salida.
2.  **Compilar**: Asegurarse de compilar el proyecto en el nodo de envío (`./compileDTsim.sh`).
3.  **Generar**: Ejecutar `python3 submitJobs.py`. Esto creará:
    *   `exec_macros/`: Directorio con los macros individuales para cada job.
    *   `logs/`: Directorio para logs de salida y error (organizados por dataset).
    *   `run_wrapper.sh`: Script wrapper para ejecutar en los nodos.
    *   `submit.sub`: Archivo de envío de HTCondor.
4.  **Enviar**: Ejecutar `condor_submit submit.sub`.

```bash
./compileDTsim.sh        # Compilar primero
python3 submitJobs.py    # Generar configuración
condor_submit submit.sub # Enviar al cluster
```

## Estructura del Proyecto

```
DTG4Simulator/
├── exampleDTSim.cc          # Main
├── include/src/             # Headers e implementaciones
├── geometry/                # Scripts Python + archivos .tg
├── macros/                  # Macros organizados (batch, interactive, settings, visualization)
├── docs/                    # Documentación técnica detallada
└── build/                   # Compilación
```

## Licencia

Geant4 BSD License - Ver [LICENSE](LICENSE)
