# Geometría del Detector

## Arquitectura: Geometría Generada desde Python

A diferencia de implementaciones tradicionales, este proyecto **NO construye la geometría directamente en C++**. En su lugar, utiliza:

1. **Script Python**: `geometry/generate_geometry_data.py`
2. **Librería mplDTs**: Extrae geometría de las cámaras DT del CMS
3. **Formato Geant4 Text Geometry**: Archivos `.tg` con sintaxis ASCII
4. **Lectura en runtime**: `G4tgbVolumeMgr` carga la geometría desde archivos

### Flujo de Generación

```
mplDTs (geometría CMS)
    ↓
generate_geometry_data.py
    ↓
geometry/geometry_concentrator.tg  ← Archivo principal
    ↓
geometry/stations/station_W-1_Sec1_St1.tg  ← Archivos individuales
geometry/stations/station_W-1_Sec1_St2.tg
...
    ↓
G4tgbVolumeMgr::ReadAndConstructDetector()  ← Carga en Geant4
```

## Script de Generación: `generate_geometry_data.py`

### Componentes Principales

El script define funciones para crear:

1. **Volúmenes de Estación** (`create_station_volume`)
2. **Yugo de Hierro** (`create_yoke`) - opcional con `include_yoke=True`
3. **Celdas de Drift** (`create_superlayer_cells`)
4. **Honeycomb de Aluminio** (`create_honeycomb`) - opcional con `include_honeycomb=True`

### Ejemplo de Uso

```python
stations_to_generate = [
    (-1, 1, 1),  # Wheel=-1, Sector=1, Station=1 (MB1)
    (-1, 1, 2),  # MB2
    (-1, 1, 3),  # MB3
    (-1, 1, 4),  # MB4
]

generate_station_geometry_ascii(
    stations_list=stations_to_generate,
    output_dir='stations',
    concentrator_template='geometry_concentrator_template',
    concentrator_output='geometry_concentrator.tg',
    include_yoke=True,
    include_honeycomb=True
)
```

**Salida**:
- `geometry/stations/station_W-1_Sec1_St{1 2 3 4}.tg` (un archivo por estación)
- `geometry/geometry_concentrator.tg` (incluye todos con `#include`)

## Formato Geant4 Text Geometry (.tg)

Los archivos `.tg` definen la geometría usando sintaxis ASCII de Geant4. Los principales elementos son:

- **Materiales**: `GasMixture` (85% Ar + 15% CO₂), `G4_Fe` (yoke), `G4_Al` (honeycomb)
- **Volúmenes**: Definidos con `:VOLU`, colocados con `:PLACE`, rotados con `:ROTM`
- **World**: Volumen madre `G4_Galactic` que contiene todas las estaciones y yokes

### Jerarquía de Volúmenes y Codificación de Celdas

Cada estación es independiente dentro del world:

```
world (G4_Galactic)
├── Station_W{wheel}_Sec{sector}_St{station} (G4_AIR)
│   ├── DriftCell_W{w}_Sec{s}_St{st}_SL{sl}_{encoding} (GasMixture)
│   │   └── [múltiples copias]
│   └── Honeycomb_W{w}_Sec{s}_St{st} (G4_Al)
└── Yoke_W{wheel}_Sec{sector}_St{station} (G4_Fe)
```

**Encoding de celdas**: El nombre del volumen incluye un código que especifica el número de celdas por capa. Por ejemplo, `DriftCell_W-1_Sec1_St2_SL1_60605959` indica 60 celdas en layer 1, 60 en layer 2, 59 en layer 3 y 59 en layer 4 (total: 238 copias del mismo volumen lógico). Este encoding es decodificado en C++ por `DriftCellSD::DecodeCellID()` para convertir el número de copia en la identificación completa (wheel, sector, station, superlayer, layer, wire).

## Regenerar Geometría

Para modificar la geometría:

```bash
cd geometry/
python generate_geometry_data.py
# Editar lista de estaciones en el script según necesidad
```

**Output**:
```
Generating DT geometry files...
Stations to generate: 12
  ✓ Generated: station_W-1_Sec1_St1.tg
  ✓ Generated: station_W-1_Sec1_St2.tg
  ...
✓ Concentrator file: geometry/geometry_concentrator.tg
  Stations: 12
  Total cells: 2847
```

Luego, no hace falta recompilar el proyecto, solo copiar los archivos `.tg` al directorio de build:

```bash
cd ../build
cmake ..
```