#include "MagneticField.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "DTSimConstants.hh"
//#include <cmath>

namespace DTSim
{

void MagneticField::GetFieldValue(const G4double point[4], G4double* field) const
{
    // Calculate radial distance in XY plane
    G4double x = point[0];
    G4double y = point[1];
    G4double r = std::sqrt(x*x + y*y);
    
    // Get field magnitude based on radius
    G4double fieldMagnitude = GetFieldMagnitude(r);
    
    // Set field components
    field[0] = 0.0;               // Bx
    field[1] = 0.0;               // By
    field[2] = fieldMagnitude;    // Bz
}

G4double MagneticField::GetFieldMagnitude(G4double r) const
{
    // Define field regions based on radial distance
    if (r < DTSim::kSolenoidRadius) {
        // Inner solenoid region - strongest field
        return DTSim::kInnerMagneticField;
    }
    else {
        // Outer detector region - medium field
        return DTSim::kOutMagneticField;
    }
}

}