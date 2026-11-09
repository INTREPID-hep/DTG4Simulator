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

La generación tiene dos mitades: ```geometry/build_geometry.py``` escribe los
```.tg``` y el concentrador, y ```geometry/detectors/<nombre>.py``` describe UN
detector. Un detector es un solo archivo y contesta tres preguntas:

    chambers()             [(etiqueta, clave, StationSpec)]: un .tg por entrada
    families()             familias de volúmenes que el concentrador pinta
    yoke(chambers, path)   archivo extra opcional (el yugo de CMS)

```build_geometry.py --detector cms``` genera las 36 cámaras de muestra, su
honeycomb y el yugo; ```--detector mutel``` genera el prototipo MuTel (dos pares
de supercapas y el muro de plomo). Para añadir un detector basta con escribir
```detectors/mi_detector.py``` y añadirlo a ```detectors/__init__.py```: no hay
que tocar el generador ni el concentrador.

El concentrador se genera a partir de las familias del propio detector, así que
un detector sin yugo ni honeycomb no necesita plantilla.

### Componentes Generados

1. **Volúmenes de Estación** - contenedor de aire para cada MB
2. **Celdas de Drift** - volumen lógico único por superlayer, múltiples copias
3. **Honeycomb de Aluminio** - capa entre superlayers, derivada del XML de CMS
4. **Yugo de Hierro** - estructura unificada con huecos para estaciones (solo CMS)

### Uso Básico

```bash
cd geometry/

# CMS: 36 cámaras de muestra + yugo + honeycomb
python build_geometry.py

# MuTel: dos pares de supercapas + muro de plomo de 5 cm
python build_geometry.py --detector mutel --output mutel     --concentrator mutel_concentrator.tg
```

El MuTel se simula apuntando el macro a su concentrador, antes de
```/run/initialize```:

```bash
/DTSim/detector/setGeometryFile geometry/mutel_concentrator.tg
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
python generate_geometry_data.py                 # = build_geometry.py --detector cms
```

La geometría se regenera sin necesidad de recompilar el proyecto. Los archivos
```.tg``` se copian al directorio de build:

```bash
cp -r stations yoke.tg geometry_concentrator.tg ../build/geometry/
```
