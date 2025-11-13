#ifndef DTSimDetectorConstruction_h
#define DTSimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4tgrMessenger;

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
    std::vector<G4LogicalVolume*> driftCellsLogicals;

    void SetVisualizationAttributes();
    
    G4tgrMessenger* fMessenger;
};

}

#endif
