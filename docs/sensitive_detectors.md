# Sensitive Detectors y Hits

## Arquitectura de Detección

El proyecto usa un único Sensitive Detector (`DriftCellSD`) para todas las celdas de drift, asignado a todos los volúmenes lógicos cuyo nombre contiene `"DriftCell"`. Los hits se recolectan en una sola colección por evento (`DriftCellHitsCollection`).

## Clase DriftCellHit

Cada hit almacena: evento ID, PDG de la partícula, carga, identificación de la celda (CellID), posición local, tiempo de drift, y energía depositada. El struct `CellID` contiene la jerarquía completa: wheel (-2 a +2), sector (1-14), station (1-4), superlayer (1-3), layer (1-4), y wire (~1-60). Incluye validación de rangos mediante `isValid()`.

## Procesamiento de Hits: DriftCellSD::ProcessHits()

El método `ProcessHits()` ejecuta para cada step en volúmenes sensibles:

1. **Filtro**: Solo registra steps con deposición de energía (edep > 0)
2. **Información geométrica**: Extrae nombre del volumen, número de copia, y posiciones (global y local)
3. **Decodificación**: Convierte nombre del volumen y copyNo en CellID completo
4. **Drift time**: Calcula tiempo incluyendo corrección de deriva (t = t₀ + |x|/v_drift, con v_drift = 54 μm/ns)
5. **Registro**: Crea objeto DriftCellHit y lo inserta en la colección

### Decodificación del CellID

El nombre del volumen tiene formato `DriftCell_W{wheel}_Sec{sector}_St{station}_SL{sl}_{encoding}`, donde el encoding especifica celdas por capa (ej: "60605959" = 60, 60, 59, 59 celdas en 4 capas). El copyNo se mapea secuencialmente a través de las capas para determinar layer y wire. Ejemplo: encoding [60, 60, 59, 59], copyNo=125 → layer=3, wire=5.

## Salida de Datos

### Formato de Hits

Cada hit imprime: evento, PDG, carga, CellID completo, posición local, tiempo de drift, y energía depositada. Ejemplo: `DriftCellHit: Event 42 PDG=13 q=-1 CellID=W-1_Sec1_St2_SL1_L3_Wire25 LocalPos=(1.2,5.3,0.1) TimeDrift=350 ns Edep=2.5 keV`

Al final del evento, los hits se guardan en el NTuple ROOT (ver [analysis.md](analysis.md) para detalles de la estructura de datos).

## Visualización de Hits

En modo interactivo con OpenGL, usar `/vis/scene/add/hits` para dibujar círculos rojos en las posiciones de los hits.

## Múltiples Hits por Celda

El sistema registra múltiples hits en la misma celda si: la partícula atraviesa en múltiples steps, diferentes partículas del shower golpean la celda, o una partícula entra/sale/reingresa.

## Validación del CellID

El método `CellID::isValid()` verifica que todos los campos estén dentro de los rangos permitidos (wheel: -2 a 2, sector: 1-14, station: 1-4, superlayer: 1-3, layer: 1-4, wire ≥ 0). Genera warnings si la decodificación falla.
