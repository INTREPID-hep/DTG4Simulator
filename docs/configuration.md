# Guía de Configuración Runtime

Esta guía documenta todos los comandos de UI (messenger) disponibles para configurar la simulación sin recompilar el código.

## Estructura de Comandos

Todos los comandos del simulador están organizados bajo el prefijo `/DTSim/`:

```
/DTSim/
├── detector/          # Configuración del detector (antes de /run/initialize)
│   ├── setGeometryFile
│   ├── enableDriftSD
│   ├── enableStationSD
│   ├── useBField
│   └── BField/
│       ├── setGlobal
│       └── setYoke
├── cellSD/            # Parámetros físicos del detector (después de /run/initialize)
│   ├── setDriftVelocity
│   ├── setMinEnergy
│   ├── setBarrierEnergy
│   ├── setWallLoss
│   ├── enableElectrostaticConfinement
│   └── enableWallCrossing
├── digitizer/         # Parámetros de digitalización (después de /run/initialize)
│   ├── setEfficiency
│   ├── setTimeResolution
│   └── setTDCResolution
└── generator/         # Configuración del generador primario
    ├── randomizePrimary
    ├── momentum
    ├── sigmaMomentum
    └── sigmaAngle
```

---

## Configuración del Detector

**IMPORTANTE**: Estos comandos deben ejecutarse **ANTES** de `/run/initialize`.

### Geometría

#### `/DTSim/detector/setGeometryFile <path>`
Especifica el archivo de geometría Text Geometry (.tg) a cargar.

**Tipo**: String  
**Valor por defecto**: `geometry/geometry_concentrator.tg`  
**Ejemplo**:
```bash
/DTSim/detector/setGeometryFile geometry/alternative_geometry.tg
```

### Sensitive Detectors

#### `/DTSim/detector/enableDriftSD <bool>`
Activa/desactiva el Sensitive Detector de las celdas de drift.

**Tipo**: Boolean  
**Valor por defecto**: `true`  
**Uso**: Desactivar para simulaciones puramente geométricas o para reducir salida de datos.

**Ejemplo**:
```bash
/DTSim/detector/enableDriftSD false
```

#### `/DTSim/detector/enableStationSD <bool>`
Activa/desactiva el Sensitive Detector de estaciones (truth segments de muones).

**Tipo**: Boolean  
**Valor por defecto**: `true`  
**Ejemplo**:
```bash
/DTSim/detector/enableStationSD false
```

### Campo Magnético

#### `/DTSim/detector/useBField <bool>`
Activa/desactiva el campo magnético en la simulación.

**Tipo**: Boolean  
**Valor por defecto**: `true`  
**Ejemplo**:
```bash
/DTSim/detector/useBField true
```

#### `/DTSim/detector/BField/setGlobal <Bx> <By> <Bz> <unit>`
Establece el vector de campo magnético global (aplicado a todo el World).

**Tipo**: G4ThreeVector con unidades  
**Valor por defecto**: `(0, 0, 0) tesla`  
**Ejemplo**:
```bash
/DTSim/detector/BField/setGlobal 0 0 3.8 tesla    # Solenoide CMS
```

#### `/DTSim/detector/BField/setYoke <Bx> <By> <Bz> <unit>`
Establece el vector de campo magnético en el yoke de hierro.

**Tipo**: G4ThreeVector con unidades  
**Valor por defecto**: `(0, 0, -2.0) tesla`  
**Ejemplo**:
```bash
/DTSim/detector/BField/setYoke 0 0 -1.8 tesla    # Campo de retorno reducido
```

---

## Parámetros Físicos del Detector

**IMPORTANTE**: Estos comandos deben ejecutarse **DESPUÉS** de `/run/initialize`.

### `/DTSim/cellSD/setDriftVelocity <value> <unit>`
Establece la velocidad de deriva de los electrones en el gas.

**Tipo**: Double con unidades de velocidad  
**Valor por defecto**: `0.054 mm/ns` (54 μm/ns para Ar-CO2 85:15)  
**Rango típico**: 40-60 μm/ns dependiendo de la mezcla de gas  
**Ejemplo**:
```bash
/DTSim/cellSD/setDriftVelocity 0.050 mm/ns    # Gas diferente
```

### `/DTSim/cellSD/setMinEnergy <value> <unit>`
Establece el umbral mínimo de energía depositada para crear un hit.

**Tipo**: Double con unidades de energía  
**Valor por defecto**: `26.6 eV` (potencial de ionización Ar-CO2)  
**Uso**: Simula el umbral de detección del gas  
**Ejemplo**:
```bash
/DTSim/cellSD/setMinEnergy 30.0 eV    # Umbral más alto
```

### `/DTSim/cellSD/setBarrierEnergy <value> <unit>`
Establece la barrera de energía del confinamiento electrostático.

**Tipo**: Double con unidades de energía  
**Valor por defecto**: `2.1 keV`  
**Uso**: Electrones con energía menor que este valor quedan atrapados por el potencial del ánodo  
**Ejemplo**:
```bash
/DTSim/cellSD/setBarrierEnergy 2.5 keV    # Confinamiento más fuerte
```

### `/DTSim/cellSD/setWallLoss <value> <unit>`
Establece la pérdida de energía al cruzar las paredes virtuales entre celdas.

**Tipo**: Double con unidades de energía  
**Valor por defecto**: `1.0 keV`  
**Uso**: Simula la atenuación por las paredes de aluminio  
**Ejemplo**:
```bash
/DTSim/cellSD/setWallLoss 0.5 keV    # Paredes más delgadas
```

### `/DTSim/cellSD/enableElectrostaticConfinement <bool>`
Activa/desactiva el modelo de confinamiento electrostático de electrones de baja energía.

**Tipo**: Boolean  
**Valor por defecto**: `true`  
**Uso**: Cuando está activo, los electrones con energía cinética menor que `fCellBarrierEnergy` quedan atrapados en la celda y no pueden escapar al ánodo vecino. Desactivar para estudiar el comportamiento sin este efecto físico.  
**Ejemplo**:
```bash
/DTSim/cellSD/enableElectrostaticConfinement false    # Desactivar confinamiento
```

### `/DTSim/cellSD/enableWallCrossing <bool>`
Activa/desactiva el modelo de pérdida de energía en paredes virtuales.

**Tipo**: Boolean  
**Valor por defecto**: `true`  
**Uso**: Cuando está activo, los electrones y positrones pierden energía (`fWallEnergyLoss`) al cruzar los límites entre celdas. Desactivar para ignorar el efecto de las paredes de aluminio.  
**Ejemplo**:
```bash
/DTSim/cellSD/enableWallCrossing false    # Desactivar pérdida en paredes
```

---

## Parámetros de Digitalización

**IMPORTANTE**: Estos comandos deben ejecutarse **DESPUÉS** de `/run/initialize`.

### `/DTSim/digitizer/setEfficiency <value>`
Establece la eficiencia de detección (probabilidad de que un hit se digitalice).

**Tipo**: Double  
**Rango**: 0.0 a 1.0  
**Valor por defecto**: `1.0` (100% de eficiencia)  
**Ejemplo**:
```bash
/DTSim/digitizer/setEfficiency 0.95    # 95% de eficiencia
```

### `/DTSim/digitizer/setTimeResolution <value> <unit>`
Establece la resolución temporal (sigma del smearing gaussiano).

**Tipo**: Double con unidades de tiempo  
**Valor por defecto**: `2.0 ns`  
**Uso**: Simula la incertidumbre temporal del detector  
**Ejemplo**:
```bash
/DTSim/digitizer/setTimeResolution 1.5 ns    # Mejor resolución
```

### `/DTSim/digitizer/setTDCResolution <value> <unit>`
Establece la resolución del TDC (tamaño del bin temporal).

**Tipo**: Double con unidades de tiempo  
**Valor por defecto**: `0.78125 ns` (25 ns / 32)  
**Uso**: Define la granularidad de la digitalización temporal  
**Ejemplo**:
```bash
/DTSim/digitizer/setTDCResolution 1.0 ns    # Bins más grandes
```

---

## Configuración del Generador

Documentado en detalle en [primary_generation.md](primary_generation.md).

### `/DTSim/generator/randomizePrimary <bool>`
Activa/desactiva la randomización de partículas primarias.

**Valor por defecto**: `false`

### `/DTSim/generator/momentum <value> <unit>`
Establece el momento del particle gun.

**Valor por defecto**: `1000.0 GeV`

### `/DTSim/generator/sigmaMomentum <value> <unit>`
Establece la dispersión gaussiana del momento.

**Valor por defecto**: `50.0 GeV`

### `/DTSim/generator/sigmaAngle <value> <unit>`
Establece la dispersión angular.

**Valor por defecto**: `2.0 deg`

---

## Orden de Ejecución

Es crucial seguir el orden correcto al usar estos comandos en macros:

```bash
# 1. Configuración del detector (ANTES de /run/initialize)
/control/execute macros/settings/verbosity.mac
/control/execute macros/settings/detector.mac

# 2. Inicialización
/run/initialize

# 3. Parámetros de física y digitalización (DESPUÉS de /run/initialize)
/control/execute macros/settings/physics.mac

# 4. Configuración del generador y ejecución
/control/execute macros/settings/primary.mac
/run/beamOn 1000
```

---

## Archivos de Configuración Predefinidos

El proyecto incluye archivos de macro organizados en `macros/settings/`:

- **detector.mac**: Configuración de geometría, SDs y campo magnético
- **physics.mac**: Parámetros físicos del detector y digitalización
- **primary.mac**: Configuración del generador de partículas
- **verbosity.mac**: Niveles de salida y debug

Ver ejemplos completos en `macros/batch/` y `macros/interactive/`.

---

## Ejemplos de Uso

### Simulación sin campo magnético
```bash
# En detector.mac
/DTSim/detector/useBField false
```

### Detector de alta eficiencia
```bash
# En physics.mac
/DTSim/digitizer/setEfficiency 1.0
/DTSim/digitizer/setTimeResolution 1.0 ns
```

### Umbral de energía más alto
```bash
# En physics.mac
/DTSim/cellSD/setMinEnergy 50.0 eV
```

### Solo hits, sin digis
```bash
# Comentar en EventAction::EndOfEvent o desactivar digitizer
# O modificar la eficiencia a 0
/DTSim/digitizer/setEfficiency 0.0
```

---

## Verificación de Configuración

Para verificar que los comandos están disponibles después de inicializar:

```bash
# En un macro o en modo interactivo
/run/initialize
/help /DTSim/
```

Esto mostrará todos los comandos disponibles bajo el directorio `/DTSim/`.
