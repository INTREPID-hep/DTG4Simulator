# Physics List

## Lista Utilizada: FTFP_BERT

El proyecto utiliza `FTFP_BERT`, la physics list por defecto de Geant4, recomendada para física de colisionadores y aplicaciones con rayos cósmicos. Ver [documentación oficial](https://geant4-userdoc.web.cern.ch/UsersGuides/PhysicsListGuide/html/reference_PL/FTFP_BERT.html).

```cpp
// En exampleDTSim.cc
auto physicsList = new FTFP_BERT;
physicsList->RegisterPhysics(new G4StepLimiterPhysics());
runManager->SetUserInitialization(physicsList);
```

### Componentes Físicos

FTFP_BERT incluye procesos hadrónicos (elásticos, inelásticos y captura) para protones, neutrones, piones, kaones, hiperones y antipartículas, combinando el modelo de Bertini (0-6 GeV) con el modelo de Fritiof (3-100 TeV). Los procesos electromagnéticos standard cubren ionización, bremsstrahlung, scattering múltiple, producción de pares y aniquilación para fotones, leptones cargados y hadrones/iones. También incluye decaimientos de partículas de vida larga.

Adicionalmente, se registra `G4StepLimiterPhysics` para controlar el tamaño de paso del tracking (?).

**No incluye**: Fotones ópticos, física de muy alta energía (>10 TeV), ni modelos de precisión para neutrones térmicos (disponibles en variantes como FTFP_BERT_HP).

## Campo Magnético

### Configuración Runtime

El campo magnético se configura mediante comandos UI en macros, sin necesidad de recompilación:

```bash
# En macros/settings/detector.mac
/DTSim/detector/useBField true                    # Habilitar campo magnético
/DTSim/detector/BField/setGlobal 0 0 0 tesla      # Campo global (x, y, z)
/DTSim/detector/BField/setYoke 0 0 -2.0 tesla     # Campo en yoke (x, y, z)
```

**Valores por defecto** (de `DTSimConstants.hh`):
- Global: (0, 0, 0) T - sin campo fuera del yoke
- Yoke: (0, 0, -2.0) T - retorno del flujo del solenoide CMS en dirección -Z

### Arquitectura del Campo

El proyecto implementa campos uniformes en dos regiones usando `G4FieldBuilder`:

1. **Campo Global**: Aplicado a todo el volumen World (típicamente cero)
2. **Campo del Yoke**: Aplicado solo a volúmenes de hierro con nombre "Yoke*"

Los comandos estándar de Geant4 `/field/` también están disponibles para ajustes avanzados.

### Ejemplo: Desactivar Campo

```bash
/DTSim/detector/useBField false
```