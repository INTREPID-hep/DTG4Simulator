#ifndef DTSimDetectorConstruction_h
#define DTSimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GenericMessenger;

namespace DTSim
{

/// Detector construction

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    void DefineCommands();
    
    std::vector<G4LogicalVolume*> fDriftCellsLogicals;
    std::vector<G4LogicalVolume*> fYokeLogicals;
    std::vector<G4LogicalVolume*> fStationLogicals;
    G4GenericMessenger* fDetMessenger;
    
    // Configuration parameters
    G4bool fUseBField;
    G4ThreeVector fGlobalField;
    G4ThreeVector fYokeField;
    G4String fGeometryFileName;
    G4bool fEnableDriftSD;
    G4bool fEnableStationSD;
};

}

#endif
