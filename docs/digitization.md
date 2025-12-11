# Digitalización (Digitization)

La digitalización es el proceso de convertir los hits de Geant4 (energía depositada, tiempo exacto) en señales digitales simuladas (Digis) que imitan la respuesta de la electrónica de lectura del detector.

## Módulo: DriftCellDigitizer

La clase `DriftCellDigitizer` (hereda de `G4VDigitizerModule`) es responsable de este proceso. Se ejecuta al final de cada evento, tomando como entrada la colección de hits (`DriftCellHitsCollection`) y produciendo una colección de digis (`DriftCellDigiCollection`).

### Proceso de Digitalización

El algoritmo de digitalización agrupa los hits por celda y selecciona solo el primero en llegar (leading edge):

1.  **Agrupamiento**: Todos los hits de Geant4 en la misma celda (mismo `CellID`) se agrupan.
2.  **Selección Temporal**: Se ordenan por tiempo de deriva (`driftTime`) y se selecciona el hit más temprano. Esto simula la electrónica que dispara con la primera ionización que llega al hilo.
3.  **Eficiencia de Detección**: Se aplica una probabilidad de detección (`kEfficiency`).
4.  **Resolución Temporal (Smearing)**:
    *   $t_{smeared} = t_{drift} + \mathcal{N}(0, \sigma_{time})$
5.  **Conversión a TDC**:
    *   $TDC = \text{int}(t_{smeared} / \text{TDC\_Resolution})$

### Parámetros de Configuración

Los parámetros de digitalización se pueden ajustar mediante comandos UI (DESPUÉS de `/run/initialize`):

```bash
# En macros/settings/physics.mac
/DTSim/digitizer/setEfficiency 1.0            # Eficiencia de detección (0.0-1.0)
/DTSim/digitizer/setTimeResolution 2.0 ns     # Sigma gaussiano del smearing
/DTSim/digitizer/setTDCResolution 0.78125 ns  # Resolución del TDC (25ns/32)
```

**Valores por defecto** (de `DTSimConstants.hh`):
- Eficiencia: 1.0 (100% de detección)
- Resolución temporal: 2.0 ns (sigma del smearing gaussiano)
- Resolución TDC: 0.78125 ns = 25ns/32 (bin estándar del TDC del CMS)

## Salida: DriftCellDigi

Cada objeto `DriftCellDigi` representa una señal digitalizada y contiene:

*   **Event ID**: Identificador del evento.
*   **Cell ID**: Identificador único de la celda (Wheel, Sector, Station, SL, Layer, Wire).
*   **TDC**: Valor digital del tiempo.
*   **Global Position**: Posición central del wire (útil para visualización, aunque en la realidad electrónica no se tiene).

### Almacenamiento en NTuple

Los digis se almacenan en el árbol `DTG4Tree` bajo el prefijo `digi_`. Ver [docs/analysis.md](analysis.md) para la lista completa de ramas.
