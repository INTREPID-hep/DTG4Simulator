# Primary Generator Action

## Configuración del Particle Gun

El `PrimaryGeneratorAction` soporta 5 tipos de partículas: e⁺, μ⁺, π⁺, K⁺, y protón. La partícula por defecto es μ⁺ (muón positivo). La selección se realiza mediante comandos UI (`/gun/particle`) o modo aleatorio.

## Parámetros del Haz

Valores por defecto: momento nominal 1000 GeV con dispersión de 50 GeV, ángulo polar θ = 90° (plano XY), ángulo azimutal φ = 0° (+X), divergencia angular de 2° en ambas direcciones, posición de origen en (0, 0, -2.5 m) sin dispersión espacial.
## Cinemática

El momento se muestrea con distribución uniforme en [p - σ_p/2, p + σ_p/2]. La energía cinética se calcula usando la relación relativista E² = p² + m², donde m es la masa de la partícula.

## Dirección del Haz

La dirección del haz se controla mediante coordenadas esféricas con dos ángulos independientes:

- **θ (theta)**: Ángulo polar desde el eje +Z (0° = +Z, 90° = plano XY, 180° = -Z)
- **φ (phi)**: Ángulo azimutal alrededor del eje Z (0° = +X, 90° = +Y, 180° = -X, 270° = -Y)

Cada ángulo tiene su propia dispersión configurable (σ_θ y σ_φ), muestreada uniformemente en [θ - σ_θ/2, θ + σ_θ/2] y [φ - σ_φ/2, φ + σ_φ/2]. Esto permite control completo de la dirección del haz en 3D y su divergencia.

El vector de dirección se calcula como:
- **d** = (sin(θ)cos(φ), sin(θ)sin(φ), cos(θ))

## Posición del Haz

La posición inicial de las partículas se puede configurar mediante un vector 3D con dispersión gaussiana opcional en cada coordenada. La posición real se calcula como: **r** = **r₀** + **N**(0, **σ**), donde **r₀** es la posición nominal y **σ** es el vector de dispersiones en (x, y, z). Esto permite simular haces con extensión espacial realista.

## Modo Aleatorio

Cuando `fRandomizePrimary = true` (default), el generador selecciona aleatoriamente entre las 5 partículas disponibles con probabilidad uniforme (20% cada una). Si es `false`, usa la partícula especificada con `/gun/particle`.

## Comandos UI

Comandos disponibles para configurar el generador primario:

- `/DTSim/generator/momentum [value] [unit]`: Momento nominal (GeV, MeV, keV)
- `/DTSim/generator/sigmaMomentum [value] [unit]`: Dispersión del momento (MeV, GeV)
- `/DTSim/generator/theta [value] [unit]`: Ángulo polar (0-180 deg)
- `/DTSim/generator/phi [value] [unit]`: Ángulo azimutal (-180 a 360 deg)
- `/DTSim/generator/sigmaTheta [value] [unit]`: Dispersión del ángulo polar (deg, rad, mrad)
- `/DTSim/generator/sigmaPhi [value] [unit]`: Dispersión del ángulo azimutal (deg, rad, mrad)
- `/DTSim/generator/position [x] [y] [z] [unit]`: Posición inicial del haz (m, cm, mm)
- `/DTSim/generator/sigmaPosition [sx] [sy] [sz] [unit]`: Dispersión espacial gaussiana (cm, mm)
- `/DTSim/generator/randomizePrimary [true|false]`: Activar/desactivar selección aleatoria de partículas

Ejemplos de uso: 
```bash
/DTSim/generator/momentum 200 GeV
/DTSim/generator/theta 90 deg      # Plano XY (perpendicular a Z)
/DTSim/generator/phi 0 deg         # Dirección hacia +X
/DTSim/generator/sigmaTheta 2 deg
/DTSim/generator/sigmaPhi 2 deg
/DTSim/generator/position 0 0 -0.25 m
/DTSim/generator/sigmaPosition 1.0 1.0 0.0 cm
/DTSim/generator/randomizePrimary false
/gun/particle mu-

# Ejemplo: haz a 45° inclinado hacia +X
/DTSim/generator/theta 45 deg
/DTSim/generator/phi 0 deg
```