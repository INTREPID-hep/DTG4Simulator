# Digitalización (Digitization)

La digitalización es el proceso de convertir los hits de Geant4 (energía depositada, tiempo exacto) en señales digitales simuladas (Digis) que imitan la respuesta de la electrónica de lectura del detector.

## Módulo: DriftCellDigitizer

La clase `DriftCellDigitizer` (hereda de `G4VDigitizerModule`) es responsable de este proceso. Se ejecuta al final de cada evento, tomando como entrada la colección de hits (`DriftCellHitsCollection`) y produciendo una colección de digis (`DriftCellDigiCollection`).

### Proceso de Digitalización

El algoritmo de digitalización aplica los siguientes pasos a cada hit:

1.  **Eficiencia de Detección**:
    *   Se aplica una probabilidad de detección (`kEfficiency`).
    *   Si el hit no pasa este corte, se descarta (simulando ineficiencia de la celda).

2.  **Resolución Temporal (Smearing)**:
    *   Se toma el tiempo de deriva (`driftTime`) del hit.
    *   Se añade un error aleatorio gaussiano (`G4RandGauss`) basado en la resolución temporal del detector (`kTimeResolution`).
    *   $t_{smeared} = t_{drift} + \mathcal{N}(0, \sigma_{time})$

3.  **Conversión a TDC**:
    *   El tiempo "manchado" (smeared) se convierte a cuentas de TDC (Time-to-Digital Converter).
    *   $TDC = \text{int}(t_{smeared} / \text{TDC\_Resolution})$
    *   Se asegura que el tiempo no sea negativo.

### Parámetros de Configuración

Los parámetros clave están definidos en `DTSimConstants.hh` (o similar, verificar implementación):

*   **Eficiencia**: Probabilidad de que un hit genere señal.
*   **Resolución Temporal**: Sigma de la gaussiana para el smearing del tiempo.
*   **Resolución TDC**: Factor de conversión de nanosegundos a cuentas TDC.

## Salida: DriftCellDigi

Cada objeto `DriftCellDigi` representa una señal digitalizada y contiene:

*   **Event ID**: Identificador del evento.
*   **Cell ID**: Identificador único de la celda (Wheel, Sector, Station, SL, Layer, Wire).
*   **TDC**: Valor digital del tiempo.
*   **Global Position**: Posición central del wire (útil para visualización, aunque en la realidad electrónica no se tiene).

### Almacenamiento en NTuple

Los digis se almacenan en el árbol `DTG4Tree` bajo el prefijo `digi_`. Ver [docs/analysis.md](analysis.md) para la lista completa de ramas.
