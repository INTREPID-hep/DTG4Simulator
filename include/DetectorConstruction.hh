#ifndef DTSimDetectorConstruction_h
#define DTSimDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
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

  private:
    void SetVisualizationAttributes();
    
    G4tgrMessenger* fMessenger;
};

}

#endif
