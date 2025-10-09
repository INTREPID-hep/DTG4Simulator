#ifndef DTSimDetectorConstruction_h
#define DTSimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;

namespace DTSim
{

/// Detector construction

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct() override;

  private:
    void ConstructMaterials();
    G4VPhysicalVolume* ConstructWorld();

    void SetVisualizationAttributes();

    // Helper method to safely get materials
    G4Material* GetMaterial(const G4String& name);

    // Store logical volumes for visualization
    G4LogicalVolume* fWorldLogical = nullptr;
    G4LogicalVolume* fDTStationLogical = nullptr;

};

}

#endif
