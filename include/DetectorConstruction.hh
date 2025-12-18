#ifndef DTSimDetectorConstruction_h
#define DTSimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4GenericMessenger;
class G4UserLimits;

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
    
    // Station User Limits
    G4double fStationStepMax;
    G4double fStationTrakMax;
    G4double fStationTimeMax;
    G4double fStationEkinMin;
    G4double fStationRangMin;
    // Yoke User Limits
    G4double fYokeStepMax;
    G4double fYokeTrakMax;
    G4double fYokeTimeMax;
    G4double fYokeEkinMin;
    G4double fYokeRangMin;

    G4UserLimits* fStationLimits;
    G4UserLimits* fYokeLimits;
};

}

#endif
