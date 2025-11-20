# Primary Generator Action

## Configuración del Particle Gun

El `PrimaryGeneratorAction` soporta 5 tipos de partículas: e⁺, μ⁺, π⁺, K⁺, y protón. La partícula por defecto es μ⁺ (muón positivo). La selección se realiza mediante comandos UI (`/gun/particle`) o modo aleatorio.

## Parámetros del Haz

Valores por defecto: momento nominal 1000 MeV con dispersión de 50 MeV, divergencia angular de 10°, posición de origen en (0, 0, -2.5 m), y dirección en el plano XY (perpendicular al eje Z del haz).
## Cinemática

El momento se muestrea con distribución uniforme en [p - σ_p/2, p + σ_p/2]. La energía cinética se calcula usando la relación relativista E² = p² + m², donde m es la masa de la partícula.

## Dirección del Haz

El ángulo se muestrea uniformemente en [-σ_angle/2, +σ_angle/2] y la dirección se establece en el plano XY. Esto simula partículas provenientes de interacciones perpendiculares al eje del haz (como en colisiones pp del CMS).

## Modo Aleatorio

Cuando `fRandomizePrimary = true` (default), el generador selecciona aleatoriamente entre las 5 partículas disponibles con probabilidad uniforme (20% cada una). Si es `false`, usa la partícula especificada con `/gun/particle`.

## Comandos UI

Comandos disponibles para configurar el generador primario:

- `/DTSim/generator/momentum [value] [unit]`: Momento nominal (GeV, MeV, keV)
- `/DTSim/generator/sigmaMomentum [value] [unit]`: Dispersión del momento (MeV, GeV)
- `/DTSim/generator/sigmaAngle [value] [unit]`: Divergencia angular (deg, rad, mrad)
- `/DTSim/generator/randomizePrimary [true|false]`: Activar/desactivar selección aleatoria de partículas

Ejemplos de uso: `/DTSim/generator/momentum 200 GeV`, `/DTSim/generator/sigmaAngle 2 deg`, `/DTSim/generator/randomizePrimary 0` seguido de `/gun/particle mu-`.