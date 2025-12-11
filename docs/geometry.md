# Geometría del Detector

## Arquitectura: Geometría Generada desde Python

Este proyecto **NO construye la geometría directamente en C++**. En su lugar utiliza:

1. **Script Python** (`geometry/generate_geometry_data.py`)
2. **Librería mplDTs** - extrae geometría CMS real de las cámaras DT
3. **Archivos Geant4 Text Geometry** (`.tg`) - formato ASCII
4. **Lectura en runtime** - `G4tgbVolumeMgr` carga la geometría

### Flujo de Generación

```
mplDTs (geometría CMS)
    ↓
generate_geometry_data.py
    ↓
Archivos .tg (uno por estación + yoke unificado)
    ↓
geometry_concentrator.tg (incluye todos)
    ↓
G4tgbVolumeMgr (carga en Geant4)
```

## Configuración Runtime

### Selección de Archivo de Geometría

El archivo de geometría se puede cambiar sin recompilar usando comandos UI:

```bash
# En macros/settings/detector.mac (ANTES de /run/initialize)
/DTSim/detector/setGeometryFile geometry/geometry_concentrator.tg
```

**Archivo por defecto**: `geometry/geometry_concentrator.tg` (definido en `DTSimConstants.hh`)

### Control de Sensitive Detectors

Los SDs se pueden activar/desactivar individualmente:

```bash
# En macros/settings/detector.mac (ANTES de /run/initialize)
/DTSim/detector/enableDriftSD true      # Hits en celdas de drift
/DTSim/detector/enableStationSD true    # Truth segments por estación
```

**Uso**: Desactivar DriftSD para simulaciones solo geométricas o para reducir salida de datos.

## Script de Generación

El script `generate_geometry_data.py` genera automáticamente:

- **Estaciones DT**: Un archivo `.tg` por estación con celdas de drift y honeycomb
- **Yoke unificado**: Archivo `yoke.tg` con yugo de hierro completo (operaciones booleanas)
- **Archivo concentrador**: Incluye todos los componentes

### Componentes Generados

1. **Volúmenes de Estación** - contenedor de aire para cada MB
2. **Celdas de Drift** - volumen lógico único por superlayer, múltiples copias
3. **Honeycomb de Aluminio** - capa entre superlayers (opcional)
4. **Yugo de Hierro** - estructura unificada con huecos para estaciones (opcional)

### Uso Básico

```python
stations_to_generate = [
    (-1, 1, 1),  # Wheel, Sector, Station (MB1)
    (-1, 1, 2),  # MB2
    # ... más estaciones
]

generate_station_geometry_ascii(
    stations_list=stations_to_generate,
    output_dir='stations',
    concentrator_template='geometry_concentrator_template',
    concentrator_output='geometry_concentrator.tg',
    include_yoke=True,        # Genera yugo de hierro
    include_honeycomb=True    # Añade honeycomb
)
```

## Archivos Generados

**Estructura de salida:**
```
geometry/
├── geometry_concentrator.tg      ← Archivo principal
├── yoke.tg                        ← Yugo unificado (TUBS + booleanas)
└── stations/
    ├── station_W-1_Sec1_St1.tg   ← Estación individual
    ├── station_W-1_Sec1_St2.tg
    └── ...
```

## Jerarquía de Volúmenes

```
world (G4_Galactic)
├── Station_W{w}_Sec{s}_St{st} (G4_AIR)
│   ├── DriftCell_..._SL{sl}_{encoding} (GasMixture)
│   └── Honeycomb_... (G4_Al)
└── Yoke_Fused (G4_Fe) ← Yugo unificado con huecos
```

**Encoding de celdas**: El nombre del volumen incluye el número de celdas por layer (ej: `60605959` = 60, 60, 59, 59 celdas). Esto permite identificar cada celda correctamente en el sensitive detector.

## Regenerar Geometría

```bash
cd geometry/
python generate_geometry_data.py
```

La geometría se regenera sin necesidad de recompilar el proyecto. Los archivos `.tg` se copian automáticamente al directorio de build.