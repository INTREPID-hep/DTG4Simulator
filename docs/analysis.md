# Análisis de Datos

## Sistema de Análisis: G4AnalysisManager

El proyecto usa `G4AnalysisManager` para salida en formato ROOT. El NTuple se almacena en el directorio `DTG4SimNTuple` y soporta merge automático de threads.

## NTuple: DTG4Tree

El árbol contiene 14 columnas con información completa de cada hit:

| ID | Nombre | Tipo | Unidad | Descripción |
|----|--------|------|--------|-------------|
| 0 | `g4dtsimHit_eventNumber` | I | - | Número de evento |
| 1 | `g4dtsimHit_PDG` | I | - | Código PDG de la partícula |
| 2 | `g4dtsimHit_q` | I | e | Carga eléctrica |
| 3 | `g4dtsimHit_wheel` | I | - | Wheel (-2 a +2) |
| 4 | `g4dtsimHit_sector` | I | - | Sector (1 a 14) |
| 5 | `g4dtsimHit_station` | I | - | Station (1 a 4) |
| 6 | `g4dtsimHit_superlayer` | I | - | SuperLayer (1 a 3) |
| 7 | `g4dtsimHit_layer` | I | - | Layer (1 a 4) |
| 8 | `g4dtsimHit_cell` | I | - | Wire (1 a ~60) |
| 9 | `g4dtsimHit_xlocal` | D | mm | Posición local X (drift) |
| 10 | `g4dtsimHit_ylocal` | D | mm | Posición local Y (wire) |
| 11 | `g4dtsimHit_zlocal` | D | mm | Posición local Z (capa) |
| 12 | `g4dtsimHit_timewithdrift` | D | ns | Tiempo con corrección de drift |
| 13 | `g4dtsimHit_edep` | D | MeV | Energía depositada |

**Tipos**:
- `I`: Integer (G4int)
- `D`: Double (G4double)

### Códigos PDG Comunes

| PDG | Partícula | Descripción |
|-----|-----------|-------------|
| 11 | e⁻ | Electrón |
| -11 | e⁺ | Positrón |
| 13 | μ⁻ | Muón negativo |
| -13 | μ⁺ | Muón positivo |
| 22 | γ | Fotón |
| 211 | π⁺ | Pión positivo |
| -211 | π⁻ | Pión negativo |
| 2212 | p | Protón |
| 2112 | n | Neutrón |

**Referencia completa**: [PDG Monte Carlo Particle Numbering Scheme](https://pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf)

## Llenado del NTuple

Cada hit en `DriftCellSD::EndOfEvent()` genera una fila en el NTuple, llenando las 14 columnas con la información del hit (evento, PDG, carga, geometría, posición, tiempo, energía).

## Archivos de Salida

Los archivos ROOT se nombran `DTG4Simulation_{runID}.root` (ej: `DTG4Simulation_0.root`). Internamente contienen el directorio `DTG4SimNTuple/` con el TTree `DTG4Tree` y sus 14 branches.