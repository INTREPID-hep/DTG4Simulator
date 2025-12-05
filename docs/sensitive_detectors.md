# Sensitive Detectors y Hits

## Arquitectura de Detección

El proyecto usa un único Sensitive Detector (`DriftCellSD`) para todas las celdas de drift, asignado a todos los volúmenes lógicos cuyo nombre contiene `"DriftCell"`. Los hits se recolectan en una sola colección por evento (`DriftCellHitsCollection`).

## Clase DriftCellHit

Cada hit almacena: evento ID, PDG de la partícula, carga, tipo de proceso (processType), identificación de la celda (CellID), posiciones (local en marco de Station y global en World), tiempo de drift, y energía depositada. 

El struct `CellID` contiene la jerarquía completa: wheel (-2 a +2), sector (1-14), station (1-4), superlayer (1-3), layer (1-4), y wire (~1-60). Incluye validación de rangos mediante `isValid()`.

## Procesamiento de Hits: DriftCellSD::ProcessHits()

El método `ProcessHits()` ejecuta para cada step en volúmenes sensibles:

1. **Filtro inicial**: Solo procesa partículas cargadas con energía depositada mayor al umbral (kMinEnergyDeposit = 26.6 eV, energía de ionización del gas Ar-CO2 85:15)

2. **Decodificación geométrica**: 
   - Extrae nombre del volumen y copyNo del **PreStepPoint** (asegura obtener el volumen sensible correcto)
   - Decodifica CellID completo del nombre del volumen
   - Valida que el CellID sea correcto antes de continuar

3. **Cálculo de posiciones**:
   - **worldPos**: Punto medio del step = (preStepPos + postStepPos) / 2
   - **cellStationPos**: Posición en marco de coordenadas de la Station (navegando al nivel 1 en la jerarquía de touchable)
   - **cellLocalPos**: Posición en marco local de la DriftCell (incluye rotaciones específicas de cada superlayer, ej: SL2 tiene rotación de 90°)

4. **Cálculo de tiempo de drift**: 
   - Tiempo medio del step: midTime = (preTime + postTime) / 2
   - Distancia radial de deriva: r = √(x² + z²) en coordenadas locales de la celda
   - Tiempo con drift: t_drift = midTime + r / v_drift (v_drift = 54 μm/ns)

5. **Información de proceso**: Obtiene el tipo de proceso físico que causó el step (electromagnetic, hadronic, etc.)

6. **Registro**: Crea objeto DriftCellHit con toda la información y lo inserta en la colección

### Decodificación del CellID

El nombre del volumen tiene formato `DriftCell_W{wheel}_Sec{sector}_St{station}_SL{sl}_{encoding}`, donde el encoding especifica celdas por capa (ej: "60605959" = 60, 60, 59, 59 celdas en 4 capas). El copyNo se mapea secuencialmente a través de las capas para determinar layer y wire. 

Ejemplo: encoding [60, 60, 59, 59], copyNo=125 → layer=3, wire=5.

## Salida de Datos

### NTuple ROOT

La información de los hits procesados por `DriftCellSD` se almacena en el árbol `DTG4Tree` bajo el prefijo `simHit_`.

El formato de salida utiliza vectores para almacenar múltiples hits por evento. Para una descripción detallada de todas las columnas disponibles, consultar la documentación de análisis: **docs/analysis.md**.

## Múltiples Hits por Celda

El sistema registra múltiples hits en la misma celda si: la partícula atraviesa en múltiples steps, diferentes partículas del shower golpean la celda, o una partícula entra/sale/reingresa. Cada step que deposite energía genera un hit independiente.

## Validación del CellID

El método `CellID::isValid()` verifica que todos los campos estén dentro de los rangos permitidos (wheel: -2 a 2, sector: 1-14, station: 1-4, superlayer: 1-3, layer: 1-4, wire ≥ 0). Si la decodificación falla, se genera un warning y el hit no se registra (retorna false).
