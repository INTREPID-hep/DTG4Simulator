# Sensitive Detectors y Hits

## Arquitectura de Detección

El proyecto utiliza dos Sensitive Detectors para capturar información de simulación:

1. **DriftCellSD**: Asignado a volúmenes `DriftCell_*`, registra hits individuales con energía depositada en cada celda de drift
2. **StationSD**: Asignado a volúmenes `Station_*`, registra segmentos de muones que atraviesan estaciones completas

Los hits se recolectan en colecciones separadas por evento: `DriftCellHitsCollection` y `DTSegmentCollection`.

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

---

## Clase DTSegment (Muon Segments)

Los **Muon Segments** capturan la trayectoria de muones que atraviesan estaciones completas. A diferencia de los hits individuales en celdas, los segmentos representan la trayectoria promedio del muón a través de toda la estación.

Cada segmento almacena:
- **StationID**: Identificación de la estación (wheel, sector, station)
- **Posición y dirección local**: En el sistema de coordenadas de la estación
- **Posición y dirección global**: En el sistema de coordenadas mundial (CMS)
- **Puntos de entrada/salida**: Para visualización solamente

El struct `StationID` contiene: wheel (-2 a +2), sector (1-14), station (1-4). Incluye validación mediante `isValid()`.

## Procesamiento de Segmentos: StationSD::ProcessHits()

El método `ProcessHits()` de `StationSD` rastrea muones (PDG = ±13) a través de las estaciones:

1. **Detección de entrada**: 
   - Al primer hit en una estación, guarda la posición de entrada
   - Utiliza un mapa `std::map<G4int, G4ThreeVector>` indexado por TrackID
   - Almacena una copia del valor

2. **Detección de salida**:
   - Identifica cuando el muón sale de la estación (volumen siguiente es "world" o "Yoke")
   - Recupera la posición de entrada del mapa

3. **Cálculo del segmento**:
   - **Posición**: Punto medio entre entrada y salida: `pos = (entry + exit) / 2`
   - **Dirección**: Vector unitario de entrada a salida: `dir = (exit - entry).unit()`
   - Transforma ambos a coordenadas locales de la estación
4. **Decodificación del StationID**:
   - Extrae wheel, sector, station del nombre del volumen: `"Station_W0_Sec1_St1"`
   - Utiliza la función utilitaria `ExtractIntAfterToken()` para parsing

5. **Registro**: Crea objeto `DTSegment` y lo inserta en `DTSegmentCollection`

El enfoque entry/exit captura la **trayectoria promedio efectiva** que mejor representa el paso del muón a través de la estación completa, en lugar de solo la dirección instantánea al entrar.

### Gestión de Memoria

- El mapa de entrada se limpia automáticamente al final de cada evento
- Los casos donde muones entran pero no salen (mueren dentro de la estación) se manejan en `EndOfEvent()`

## Utilidad: ExtractIntAfterToken()

El módulo `DTSimUtils` proporciona funciones auxiliares para parsear nombres de volúmenes:

```cpp
G4int ExtractIntAfterToken(const G4String& str, const G4String& token, size_t startPos = 0);
```

Esta función extrae valores enteros después de tokens específicos (ej: `"_W"`, `"_Sec"`, `"_St"`), reduciendo código repetitivo de parsing. Se utiliza en ambos `DriftCellSD` y `StationSD`.

## Salida de Datos

### NTuple ROOT

La información de hits y segmentos se almacena en el árbol `DTG4Tree`:
- **SimHits** (DriftCellSD): Prefijo `simHit_` - hits individuales con energía depositada
- **Muon Segments** (StationSD): Prefijo `seg_` - trayectorias de muones por estación

El formato utiliza vectores para múltiples hits/segmentos por evento. Ver **docs/analysis.md** para descripción completa de todas las columnas.
