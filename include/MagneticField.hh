#ifndef DTSimMagneticField_hh
#define DTSimMagneticField_hh 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

namespace DTSim
{

class MagneticField : public G4MagneticField
{
  public:
      MagneticField() = default;
      virtual ~MagneticField() = default;

      // Main method - returns field value at given point
      virtual void GetFieldValue(const G4double point[4], G4double* field) const override;

  private:
      // Helper method to get field magnitude based on radial position
      G4double GetFieldMagnitude(G4double radiusXY) const;
};

}

#endif