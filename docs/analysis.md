# Análisis de Datos

## Sistema de Análisis: G4AnalysisManager

El proyecto usa `G4AnalysisManager` para salida en formato ROOT. El NTuple se almacena en el directorio `DTG4SimNTuple` y soporta merge automático de threads.

## NTuple: DTG4Tree

El árbol `DTG4Tree` almacena la información **por evento**. A diferencia de versiones anteriores, se utilizan `std::vector` para almacenar múltiples hits, digis y partículas generadas dentro de una sola entrada (fila) del TTree.

### Estructura de Ramas (Branches)

El árbol se divide en tres categorías principales:

#### 1. Información del Evento
| Nombre | Tipo | Descripción |
|--------|------|-------------|
| `event_eventNumber` | Int | Número de identificación del evento |

#### 2. Generador (Truth Level)
Prefijo: `gen_`
| Nombre | Tipo | Descripción |
|--------|------|-------------|
| `gen_nGenParts` | Int | Número de partículas primarias generadas |
| `gen_pdgId` | Vector\<Int> | Código PDG de la partícula |
| `gen_charge` | Vector\<Int> | Carga eléctrica |
| `gen_pt` | Vector\<Double> | Momento transversal ($p_T$) [GeV] |
| `gen_eta` | Vector\<Double> | Pseudorapidez ($\eta$) |
| `gen_phi` | Vector\<Double> | Ángulo azimutal ($\phi$) [rad] |

#### 3. SimHits (Geant4 Hits)
Prefijo: `simHit_`
| Nombre | Tipo | Descripción |
|--------|------|-------------|
| `simHit_nSimHits` | Int | Número total de hits en el evento |
| `simHit_PDG` | Vector\<Int> | Código PDG de la partícula que causó el hit |
| `simHit_q` | Vector\<Int> | Carga de la partícula |
| `simHit_wheel` | Vector\<Int> | ID de Wheel (-2 a +2) |
| `simHit_sector` | Vector\<Int> | ID de Sector (1 a 14) |
| `simHit_station` | Vector\<Int> | ID de Station (1 a 4) |
| `simHit_superlayer` | Vector\<Int> | ID de SuperLayer (1 a 3) |
| `simHit_layer` | Vector\<Int> | ID de Layer (1 a 4) |
| `simHit_cell` | Vector\<Int> | Número de Wire |
| `simHit_xlocal` | Vector\<Double> | Posición X local (dirección de drift) [mm] |
| `simHit_ylocal` | Vector\<Double> | Posición Y local (a lo largo del wire) [mm] |
| `simHit_zlocal` | Vector\<Double> | Posición Z local [mm] |
| `simHit_time` | Vector\<Double> | Tiempo global + drift [ns] |
| `simHit_edep` | Vector\<Double> | Energía depositada [MeV] |
| `simHit_process_type` | Vector\<Int> | Tipo de proceso físico (Geant4 code) |

#### 4. Digis (Digitalización)
Prefijo: `digi_`
| Nombre | Tipo | Descripción |
|--------|------|-------------|
| `digi_nDigis` | Int | Número de digis creados |
| `digi_wheel` | Vector\<Int> | ID de Wheel |
| `digi_sector` | Vector\<Int> | ID de Sector |
| `digi_station` | Vector\<Int> | ID de Station |
| `digi_superlayer` | Vector\<Int> | ID de SuperLayer |
| `digi_layer` | Vector\<Int> | ID de Layer |
| `digi_cell` | Vector\<Int> | Número de Wire |
| `digi_TDC` | Vector\<Int> | Valor TDC (Time-to-Digital Converter) |

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

Cada evento en `RunAction::EndOfEventAction()` genera una fila en el NTuple. Los vectores se limpian al inicio de cada evento y se llenan con los hits y digis acumulados.

## Archivos de Salida

Los archivos ROOT se nombran `DTG4Simulation_{runID}.root` (ej: `DTG4Simulation_0.root`). Internamente contienen el directorio `DTG4SimNTuple/` con el TTree `DTG4Tree`.