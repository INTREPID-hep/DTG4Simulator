#ifndef DTSimConstants_hh
#define DTSimConstants_hh 1

#include "G4SystemOfUnits.hh"

namespace DTSim
{

// Drift tube physical constants
constexpr G4double kDriftVelocity = 54.0*micrometer/nanosecond;

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
constexpr G4double kInnerMagneticField = 10*tesla;
constexpr G4double kOutMagneticField = -2.0*tesla; // Negative for field direction
constexpr G4double kYokeMagneticField = -0.5*tesla;
const G4double kSolenoidRadius = 3.5*m;   // Inner solenoid radius

}

#endif