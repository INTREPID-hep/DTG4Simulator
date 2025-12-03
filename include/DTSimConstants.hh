#ifndef DTSimConstants_hh
#define DTSimConstants_hh 1

#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

namespace DTSim
{

// Drift tube physical constants
constexpr G4double kDriftVelocity = 54.0*micrometer/nanosecond;
constexpr G4double kMinEnergyDeposit = 26.6 * eV; // Minimum energy deposit to ionize a gas mixture Ar-Co2 85:15

// Geometry limits
constexpr G4int kMinWheel = -2;
constexpr G4int kMaxWheel = 2;
constexpr G4int kMinSector = 1;
constexpr G4int kMaxSector = 14;
constexpr G4int kMinStation = 1;
constexpr G4int kMaxStation = 4;
constexpr G4int kMinSuperLayer = 1;
constexpr G4int kMaxSuperLayer = 3;
constexpr G4int kMinLayer = 1;
constexpr G4int kMaxLayer = 4;

// Magnetic field constants
constexpr G4double kInnerMagneticField = 3.8 *tesla; // Inner solenoid field - currently unused
const G4double kSolenoidRadius = 3.5*m;   // Inner solenoid radius - currently unused
constexpr G4double kOutMagneticField = 0.0 *tesla; 
constexpr G4double kYokeMagneticField = -2.0 *tesla; // Negative for field direction

// Primary generator defaults
constexpr G4double kDefaultMomentum = 1000.0*MeV;
constexpr G4double kDefaultSigmaMomentum = 50.0*MeV;
constexpr G4double kDefaultSigmaAngle = 2.0*deg;
constexpr G4double kDefaultParticleEnergy = 10.0*GeV;
const G4ThreeVector kDefaultParticlePosition = G4ThreeVector(0., 0., -2.5*m);
}

#endif