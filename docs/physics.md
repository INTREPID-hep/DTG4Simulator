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

El proyecto simplifica el campo magnético como uniforme en dos regiones: el yoke de hierro con $-2.0 T$ en dirección -Z (retorno del flujo del solenoide del CMS), y el volumen externo con $0 T$ (valores modificables en DTSimConstants.hh, requiere recompilar). El campo se incluye en la simulación únicamente con el flag `-B` en línea de comandos:

```bash
./exampleDTSim -B -m run.mac
```

Los comandos UI `/field/` permiten ajustar el campo global, mientras que comandos específicos del yoke (`/field/yoke_.../`) controlan el campo local en los volúmenes de hierro.